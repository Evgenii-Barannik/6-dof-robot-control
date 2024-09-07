#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <inttypes.h>  
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/// Ensures that assertions are not disabled
#undef NDEBUG
#include <assert.h>
///

/// CONSTANTS THAT SHOULD NOT BE CHANGED
const double DEGREES_TO_RADIANS = (2 * M_PI / 360);
const int NUM_OF_HEXAGON_VERTICES = 6;
const int NUM_OF_VECTORS = NUM_OF_HEXAGON_VERTICES + 2;
const int NUM_OF_SLIDERS = 6;
///

/// CONSTANTS THAT CAN BE CHANGED
const bool PRINT_DEBUG_INFO = true;
const double NEEDLE_LENGTH = 3.0;
const double TABLE_RADIUS = 2.5; // Radius of the hexagonal table in cm.
const double TABLE_ANGLE = DEGREES_TO_RADIANS * 80; // Angle between vertices of the hexagonal table in radians. A regular hexagon will have 60 deg.
const double BASE_RADIUS = 7.0; // Radius of the hexagonal base in cm.
const double BASE_ANGLE = DEGREES_TO_RADIANS * 15; // Angle between vertices of the hexagonal base in radians. A regular hexagon will have 60 deg.
const double MAX_POSSIBLE_HEIGHT = 20.0; // Highest possible position for sliders in cm.
const double MIN_POSSIBLE_HEIGHT = 0; // Lowest possible position for sliders in cm.
const double ARM_LENGTH = 8.0; // Length of all robot arms in cm.
const double DYNAMIC_ACCURACY = 0.99;  // should be between 0.0 and 1.0 
///

// Full description of the needle position, and thus, whole system (ideally).
typedef struct { 
	double x, y, z; // Coordinates of the needle end in cm.
	double phi, theta, psi; // Euler angles for the needle in radians.
} NeedleCoordinates;

typedef struct {
    gsl_vector* coordinates[NUM_OF_VECTORS];
} Table_Model;

// Array of pointers to NULL (if slider coordinates were not found)
// or pointers to slider coordinates (if slider coordinates were found).
typedef struct {
    gsl_vector* coordinates[NUM_OF_SLIDERS];
} Sliders_Model_Nullable;

typedef struct{
	// Pointers to arrays with size equal to the number of frames  
	Table_Model* table_positions; 
	Sliders_Model_Nullable* slider_positions;
	uint64_t* timestamps;

	// Pointers to single structs
	NeedleCoordinates* initial_state;
	NeedleCoordinates* final_state;
} Trajectory;

bool all_not_NULL(const Sliders_Model_Nullable sliders) {
    for (int i = 0; i < NUM_OF_SLIDERS; ++i) {
        gsl_vector* vec_p = sliders.coordinates[i];
        if (vec_p == NULL) {
            return false;
        }
    }
    return true;
}
void print_needle_coordinates(NeedleCoordinates* table) {
	printf("Linear coordinates (x, y, z) in cm:        (%f, %f, %f)\n", table->x, table->y, table->z);
	printf("Euler angles (phi, theta, psi) in radians: (%f, %f, %f)\n", table->phi, table->theta, table->psi);
}

void print_table(const Table_Model table) {
	for(int i = 0; i < NUM_OF_VECTORS; i++) {

		printf("%d: (%f, %f, %f)", i,
				gsl_vector_get(table.coordinates[i], 0),
				gsl_vector_get(table.coordinates[i], 1),
				gsl_vector_get(table.coordinates[i], 2));

		if (i == NUM_OF_VECTORS - 2) {
			printf(" <-- table centroid");
		}
		if (i == NUM_OF_VECTORS - 1){
			printf(" <-- needle end");
		}
		printf("\n");
	}
}

void print_sliders(const Sliders_Model_Nullable sliders) {
	for(int i = 0; i < NUM_OF_SLIDERS; i++) {
		if (sliders.coordinates[i] != NULL) {
			printf("Slider %d: (%f, %f, %f)\n", i,
					gsl_vector_get(sliders.coordinates[i], 0),
					gsl_vector_get(sliders.coordinates[i], 1),
					gsl_vector_get(sliders.coordinates[i], 2));
		} else {
			printf("Slider %d: ERROR\n", i);
		}
	}
}

void print_matrix_4D(const gsl_matrix* m) {
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			printf(j==3?"%9.6f\n":"%9.6f ", gsl_matrix_get(m,i,j));
}

void print_trajectory(const Trajectory* trajectory, const size_t num_of_frames, const gsl_matrix* transformation_matrix) {
	printf("\n==== Trajectory ====");
	printf("\nNumber of frames: %d", (int) num_of_frames);
	printf("\nTrajectory start time: %" PRIu64 " microseconds\n", trajectory->timestamps[0]);
	printf("Trajectory end time:   %" PRIu64 " microseconds\n", trajectory->timestamps[num_of_frames - 1]);


	printf("\n==== Needle coordinates for initial state:\n");
	print_needle_coordinates(trajectory->initial_state);
	printf("\n==== Needle coordinates for final state:\n");
	print_needle_coordinates(trajectory->final_state);
	printf("\n==== Affine transformation matrix:\n");
	print_matrix_4D(transformation_matrix);

	for (int i = 0; i < num_of_frames; i++) {
		double progress = (double) i / (double) (num_of_frames - 1);
		printf("\n==== Frame %d / %d\nProgress: %f\nTimestamp: %" PRIu64 "\n", i, (int) num_of_frames - 1, progress, trajectory->timestamps[i]);
		printf("\nTable coordinates (x, y, z):\n");
		print_table(trajectory->table_positions[i]);
		printf("\nSlider coordinates (x, y, z):\n");
		print_sliders(trajectory->slider_positions[i]);
	}
}



gsl_matrix* get_transformation_matrix(const NeedleCoordinates state_A, const NeedleCoordinates state_B) {
	// M_0B, transformation matrix for (state_0) -> (state_B)
	gsl_matrix* M_0B = gsl_matrix_alloc(4, 4);
	gsl_matrix_set(M_0B, 0, 0, cos(state_B.psi) * cos(state_B.theta));
	gsl_matrix_set(M_0B, 0, 1, cos(state_B.psi) * sin(state_B.theta) * sin(state_B.phi) - sin(state_B.psi) * cos(state_B.phi));
	gsl_matrix_set(M_0B, 0, 2, cos(state_B.psi) * sin(state_B.theta) * cos(state_B.phi) + sin(state_B.psi) * sin(state_B.phi));
	gsl_matrix_set(M_0B, 0, 3, state_B.x);
	gsl_matrix_set(M_0B, 1, 0, sin(state_B.psi) * cos(state_B.theta));
	gsl_matrix_set(M_0B, 1, 1, sin(state_B.psi) * sin(state_B.theta) * sin(state_B.phi) + cos(state_B.psi) * cos(state_B.phi));
	gsl_matrix_set(M_0B, 1, 2, sin(state_B.psi) * sin(state_B.theta) * cos(state_B.phi) - cos(state_B.psi) * sin(state_B.phi));
	gsl_matrix_set(M_0B, 1, 3, state_B.y);
	gsl_matrix_set(M_0B, 2, 0, -sin(state_B.theta));
	gsl_matrix_set(M_0B, 2, 1, cos(state_B.theta) * sin(state_B.phi));
	gsl_matrix_set(M_0B, 2, 2, cos(state_B.theta) * cos(state_B.phi));
	gsl_matrix_set(M_0B, 2, 3, state_B.z);
	gsl_matrix_set(M_0B, 3, 0, 0.0);
	gsl_matrix_set(M_0B, 3, 1, 0.0);
	gsl_matrix_set(M_0B, 3, 2, 0.0);
	gsl_matrix_set(M_0B, 3, 3, 1.0);

	// M_0A, transformation matrix for (state_0) -> (state_A)
	gsl_matrix* M_0A = gsl_matrix_alloc(4, 4);
	gsl_matrix_set(M_0A, 0, 0, cos(state_A.psi) * cos(state_A.theta));
	gsl_matrix_set(M_0A, 0, 1, cos(state_A.psi) * sin(state_A.theta) * sin(state_A.phi) - sin(state_A.psi) * cos(state_A.phi));
	gsl_matrix_set(M_0A, 0, 2, cos(state_A.psi) * sin(state_A.theta) * cos(state_A.phi) + sin(state_A.psi) * sin(state_A.phi));
	gsl_matrix_set(M_0A, 0, 3, state_A.x);
	gsl_matrix_set(M_0A, 0, 3, state_A.x);
	gsl_matrix_set(M_0A, 1, 0, sin(state_A.psi) * cos(state_A.theta));
	gsl_matrix_set(M_0A, 1, 1, sin(state_A.psi) * sin(state_A.theta) * sin(state_A.phi) + cos(state_A.psi) * cos(state_A.phi));
	gsl_matrix_set(M_0A, 1, 2, sin(state_A.psi) * sin(state_A.theta) * cos(state_A.phi) - cos(state_A.psi) * sin(state_A.phi));
	gsl_matrix_set(M_0A, 1, 3, state_A.y);
	gsl_matrix_set(M_0A, 2, 0, -sin(state_A.theta));
	gsl_matrix_set(M_0A, 2, 1, cos(state_A.theta) * sin(state_A.phi));
	gsl_matrix_set(M_0A, 2, 2, cos(state_A.theta) * cos(state_A.phi));
	gsl_matrix_set(M_0A, 2, 3, state_A.z);
	gsl_matrix_set(M_0A, 3, 0, 0.0);
	gsl_matrix_set(M_0A, 3, 1, 0.0);
	gsl_matrix_set(M_0A, 3, 2, 0.0);
	gsl_matrix_set(M_0A, 3, 3, 1.0);

	int signum;
	gsl_permutation * p = gsl_permutation_alloc(4);
	gsl_matrix* M_A0 = gsl_matrix_alloc(4, 4);
	gsl_matrix* M_0A_temp = gsl_matrix_alloc(4, 4);
	gsl_matrix_memcpy(M_0A_temp, M_0A);
	gsl_linalg_LU_decomp (M_0A_temp, p, &signum);
	gsl_linalg_LU_invert (M_0A_temp, p, M_A0);
	gsl_permutation_free(p);

	// M_AB, transformation matrix for (state_A) -> (state_B)
	// M_AB = M_0B * (M_0A)^(-1)
	// M_AB = M_0B * M_A0
	gsl_matrix* M_AB = gsl_matrix_alloc(4, 4);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, M_0B, M_A0,  0.0, M_AB);
	return M_AB;
}

void set_transformation_frames(
	gsl_matrix** transformation_frames,
	const NeedleCoordinates state_A,
	const NeedleCoordinates state_B,
	const size_t num_of_frames
	) {
	
	const double delta_X = state_B.x - state_A.x;
	const double delta_Y = state_B.y - state_A.y;
	const double delta_Z = state_B.z - state_A.z;
	const double delta_Phi = state_B.phi - state_A.phi;
	const double delta_Theta = state_B.theta - state_A.theta;
	const double delta_Psi = state_B.psi - state_A.psi;

	assert(num_of_frames >= 2);
	assert(delta_Phi < 2*M_PI);
	assert(delta_Theta < 2*M_PI);
	assert(delta_Psi < 2*M_PI);

	for (size_t j = 0; j < num_of_frames; j++) {
		const double progress_coefficient = (double) j / (double) (num_of_frames-1);

		const double intermediate_delta_X = delta_X * progress_coefficient;
		const double intermediate_delta_Y = delta_Y * progress_coefficient;
		const double intermediate_delta_Z = delta_Z * progress_coefficient;
		const double intermediate_delta_Phi = delta_Phi * progress_coefficient;
		const double intermediate_delta_Theta = delta_Theta * progress_coefficient;
		const double intermediate_delta_Psi = delta_Psi * progress_coefficient;

		const NeedleCoordinates state_T = {
			state_A.x + intermediate_delta_X,
			state_A.y + intermediate_delta_Y,
			state_A.z + intermediate_delta_Z,
			state_A.phi + intermediate_delta_Phi,
			state_A.theta + intermediate_delta_Theta,
			state_A.psi + intermediate_delta_Psi
		};
		
		transformation_frames[j] = get_transformation_matrix(state_A, state_T); // M_AT
	}
}

// Coordinates of sliders at their highest possible position (at rest).
Sliders_Model_Nullable* get_sliders_in_state_0(const double radius, const double angle_in_radians, bool print_debug_info) {
	Sliders_Model_Nullable* sliders = (Sliders_Model_Nullable*) malloc(sizeof(Sliders_Model_Nullable));

	for (int i = 0; i < NUM_OF_SLIDERS; i++) {
		sliders->coordinates[i] = gsl_vector_alloc(3);
	}

	gsl_vector_set(sliders->coordinates[0], 0, radius * sin(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[0], 1, radius * cos(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[0], 2, MAX_POSSIBLE_HEIGHT);
	gsl_vector_set(sliders->coordinates[1], 0, radius * sin(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[1], 1, radius * cos(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[1], 2, MAX_POSSIBLE_HEIGHT);
	gsl_vector_set(sliders->coordinates[2], 0, radius * sin(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[2], 1, radius * cos(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[2], 2, MAX_POSSIBLE_HEIGHT);
	gsl_vector_set(sliders->coordinates[3], 0, radius * sin(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[3], 1, radius * cos(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[3], 2, MAX_POSSIBLE_HEIGHT);
	gsl_vector_set(sliders->coordinates[4], 0, radius * sin(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[4], 1, radius * cos(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[4], 2, MAX_POSSIBLE_HEIGHT);
	gsl_vector_set(sliders->coordinates[5], 0, radius * sin(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[5], 1, radius * cos(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[5], 2, MAX_POSSIBLE_HEIGHT);


	if (print_debug_info) {
		printf("\n==== Slider coordinates (x, y, z) in cm for state_0:\n");
		print_sliders(*sliders);
	}

	return sliders;
}

Table_Model* get_table_in_state_0(bool print_debug_info) {
	Table_Model* table = (Table_Model*)malloc(sizeof(Table_Model));
	for (int i = 0; i < NUM_OF_VECTORS; i++) {
		table->coordinates[i] = gsl_vector_alloc(3);
	}

	Sliders_Model_Nullable* sliders_in_state_0 = get_sliders_in_state_0(BASE_RADIUS, BASE_ANGLE, false);

	double radius = TABLE_RADIUS;
	double angle_in_radians = TABLE_ANGLE;
	gsl_vector_set(table->coordinates[0], 0, radius * sin(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(table->coordinates[0], 1, radius * cos(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(table->coordinates[1], 0, radius * sin(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(table->coordinates[1], 1, radius * cos(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(table->coordinates[2], 0, radius * sin(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(table->coordinates[2], 1, radius * cos(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(table->coordinates[3], 0, radius * sin(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(table->coordinates[3], 1, radius * cos(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(table->coordinates[4], 0, radius * sin(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(table->coordinates[4], 1, radius * cos(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(table->coordinates[5], 0, radius * sin(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(table->coordinates[5], 1, radius * cos(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));

	// Now we calculate the z-coordinate for the first table vertex. It will be the same for all 6 hexagon vertices and centroid.
	double slider_x = gsl_vector_get(sliders_in_state_0->coordinates[0], 0);
	double slider_y = gsl_vector_get(sliders_in_state_0->coordinates[0], 1);
	double slider_z = gsl_vector_get(sliders_in_state_0->coordinates[0], 2);
	double table_x = gsl_vector_get(table->coordinates[0], 0);
	double table_y = gsl_vector_get(table->coordinates[0], 1);
	double deltaX = slider_x - table_x;
	double deltaY = slider_y - table_y;
	double deltaZ_squared = pow(ARM_LENGTH, 2) - pow(deltaX, 2) - pow(deltaY, 2);
	assert(deltaZ_squared > 0);
	double table_z = slider_z - sqrt(deltaZ_squared); 

	for (int i = 0; i < NUM_OF_HEXAGON_VERTICES; i++) {
		gsl_vector_set(table->coordinates[i], 2, table_z);
	}

	gsl_vector_set(table->coordinates[6], 0, 0.0); 
	gsl_vector_set(table->coordinates[6], 1, 0.0);
	gsl_vector_set(table->coordinates[6], 2, table_z); 

	gsl_vector_set(table->coordinates[7], 0, 0.0);
	gsl_vector_set(table->coordinates[7], 1, 0.0);
	gsl_vector_set(table->coordinates[7], 2, table_z - NEEDLE_LENGTH);

	for (int i = 0; i < NUM_OF_SLIDERS; ++i) {
		gsl_vector_free(sliders_in_state_0->coordinates[i]);
	}
	free(sliders_in_state_0);

	if (print_debug_info) {
		printf("==== Table coordinates (x, y, z) in cm for state_0:\n");
		print_table(*table);
	}

	return table;
}

Sliders_Model_Nullable* find_sliders(const Table_Model table_transformed_vecs, const Sliders_Model_Nullable sliders_in_state_0) {
	Sliders_Model_Nullable* sliders = (Sliders_Model_Nullable*) malloc(sizeof(Sliders_Model_Nullable));

	for (int i = 0; i < NUM_OF_SLIDERS; i++) {
		double slider_x = gsl_vector_get(sliders_in_state_0.coordinates[i], 0);
		double slider_y = gsl_vector_get(sliders_in_state_0.coordinates[i], 1);

		double deltaX = gsl_vector_get(sliders_in_state_0.coordinates[i], 0) - gsl_vector_get(table_transformed_vecs.coordinates[i], 0);
		double deltaY = gsl_vector_get(sliders_in_state_0.coordinates[i], 1) - gsl_vector_get(table_transformed_vecs.coordinates[i], 1);
		double deltaZ_squared = pow(ARM_LENGTH, 2) - pow(deltaX, 2) - pow(deltaY, 2);

		if (deltaZ_squared > 0) {
			// This is one of two possible solutions for z coordinate of slider, the higher one;
			// Lower solution will have minus sqrt(targerted_deltaZ_squared) instead of plus sqrt(targerted_deltaZ_squared);
			// We don't use second solution at alle because we fear damage during possible solution switching.
			double slider_z = gsl_vector_get(table_transformed_vecs.coordinates[i], 2) + sqrt(deltaZ_squared);
			if ((slider_z < MAX_POSSIBLE_HEIGHT) && (slider_z > MIN_POSSIBLE_HEIGHT)) {
				sliders->coordinates[i] = gsl_vector_alloc(3);
				gsl_vector_set(sliders->coordinates[i], 0, slider_x);
				gsl_vector_set(sliders->coordinates[i], 1, slider_y);
				gsl_vector_set(sliders->coordinates[i], 2, slider_z);
			} else {
				sliders->coordinates[i] = NULL;
			};
		} else {
			sliders->coordinates[i] = NULL;
		}
	}
	return sliders;
}

void set_sliders_frames(
		Sliders_Model_Nullable* slider_frames,
		const size_t num_of_frames,
		const Table_Model* table_frames,
		const Sliders_Model_Nullable sliders_in_state_0
		) {

	for (size_t j = 0; j < num_of_frames; j++)  {
		Table_Model table_for_this_frame = table_frames[j];
		Sliders_Model_Nullable* sliders_for_this_frame = find_sliders(table_for_this_frame, sliders_in_state_0);
		slider_frames[j] = *sliders_for_this_frame;
	}
}

Table_Model* get_transformed_table(const gsl_matrix* T, const Table_Model initial_table){
	Table_Model* transformed_table = (Table_Model*)malloc(sizeof(Table_Model));

	for (int i = 0; i < NUM_OF_VECTORS; i++) {
		const gsl_vector* initial_vec_in_3D = initial_table.coordinates[i];

		// Uplifting to 4D
		gsl_vector* initial_vec_in_4D = gsl_vector_alloc(4);
		gsl_vector_set(initial_vec_in_4D, 0, gsl_vector_get(initial_vec_in_3D, 0)); // We take X coordinate of 3D vector
		gsl_vector_set(initial_vec_in_4D, 1, gsl_vector_get(initial_vec_in_3D, 1)); // We take Y coordinate of 3D vector
		gsl_vector_set(initial_vec_in_4D, 2, gsl_vector_get(initial_vec_in_3D, 2)); // We take Z coordinate of 3D vector
		gsl_vector_set(initial_vec_in_4D, 3, 1.0); // We use new coordinate for Tau

		// Linear transformation in 4D
		gsl_vector* transformed_vec_in_4D = gsl_vector_alloc(4);
		gsl_blas_dgemv(CblasNoTrans, 1.0, T, initial_vec_in_4D, 0.0, transformed_vec_in_4D);

		/*double point_z = gsl_vector_get(transformed_vec_in_4D, 2);*/
		/*assert(point_z >= MIN_POSSIBLE_HEIGHT);*/

		/*// Projection back to 3D*/
		transformed_table->coordinates[i] = gsl_vector_alloc(3);
		gsl_vector_set(transformed_table->coordinates[i], 0, gsl_vector_get(transformed_vec_in_4D, 0));
		gsl_vector_set(transformed_table->coordinates[i], 1, gsl_vector_get(transformed_vec_in_4D, 1));
		gsl_vector_set(transformed_table->coordinates[i], 2, gsl_vector_get(transformed_vec_in_4D, 2));

		gsl_vector_free(initial_vec_in_4D);
		gsl_vector_free(transformed_vec_in_4D);
	}
	return transformed_table;
}

void set_table_frames(Table_Model* table_frames, const size_t num_of_frames, const gsl_matrix** M_frames, const Table_Model initial_table) {
	for (size_t j = 0; j < num_of_frames; j++) {
		table_frames[j] = *get_transformed_table(M_frames[j], initial_table);
		// const gsl_matrix** prohibits this:
		// gsl_matrix_set(M_frames[j], 1, 1, 0.0);
	}
}

Trajectory* get_trajectory(const NeedleCoordinates initial_state, const NeedleCoordinates final_state, uint64_t movement_start_mks) {
	const uint64_t movement_duration_in_mks = 3 * 1000 * 1000; // TODO - calculate this
	const size_t num_of_frames = 10; // TODO - calculate this
	assert(num_of_frames >= 2);

	// Struct with array of pointers to vectors that describe Table 3D model in state_0.
	const Table_Model table_in_state_0 = *get_table_in_state_0(PRINT_DEBUG_INFO);
	double needle_z_in_state_0 = gsl_vector_get(table_in_state_0.coordinates[7], 2);
	const NeedleCoordinates state_0 = {0.0, 0.0, needle_z_in_state_0, 0.0, 0.0, 0.0};

	const gsl_matrix* M_0A = get_transformation_matrix(state_0, initial_state);

	// Struct with array of pointers to vectors that describe Table 3D model in state_A
	Table_Model table_in_state_A = *get_transformed_table(M_0A, table_in_state_0);

	// Struct with array of pointers to vectors that describe Sliders 3D model in zero state
	Sliders_Model_Nullable sliders_in_state_0 = *get_sliders_in_state_0(BASE_RADIUS, BASE_ANGLE, PRINT_DEBUG_INFO);

	// First frame is the the identity transformation (Id).
	// Last frame is the full transformation (M_AB), which brings (state_A) to the (state_B).
	// Other frames describe transformation that is "in progress", which bring (state_A) to some intermediate states (state_Tn)
	// Example for the case with 4 frames:
	// 	frame 0 is the Id;
	//	frame 1 is the M_AT1
	//	frame 2 is the M_AT2
	//	frame 3 is the M_AB;
	gsl_matrix* M_frames[num_of_frames];
	set_transformation_frames(M_frames, initial_state, final_state, num_of_frames);

	// Each 3D model is produced by applying one of transformations "in progress" to the 3D model for initial state.
	Table_Model table_frames[num_of_frames];
	set_table_frames(table_frames, num_of_frames, (const gsl_matrix**) M_frames, table_in_state_A);

	Sliders_Model_Nullable slider_frames[num_of_frames];
	set_sliders_frames(slider_frames, num_of_frames, table_frames, sliders_in_state_0);
	// print_transformed_slider_coordinates(slider_frames, num_of_frames);
	//
	uint64_t timestamps[num_of_frames];
	for (int i = 0; i < num_of_frames; i++) {
    		const double progress = (double)i / (double)(num_of_frames - 1);
	    	timestamps[i] = movement_start_mks + (uint64_t)(movement_duration_in_mks * progress);
	}

	Trajectory* trajectory = (Trajectory*)malloc(sizeof(Trajectory));
	trajectory->table_positions = (Table_Model*)malloc(num_of_frames * sizeof(Table_Model));
	trajectory->slider_positions = (Sliders_Model_Nullable*)malloc(num_of_frames * sizeof(Sliders_Model_Nullable));
	trajectory->timestamps = (uint64_t*)malloc(num_of_frames * sizeof(uint64_t));
	trajectory->initial_state = malloc(sizeof(NeedleCoordinates));
	trajectory->final_state = malloc(sizeof(NeedleCoordinates));

	memcpy(trajectory->table_positions, table_frames, num_of_frames * sizeof(Table_Model));
	memcpy(trajectory->slider_positions, slider_frames, num_of_frames * sizeof(Sliders_Model_Nullable));
	memcpy(trajectory->timestamps, timestamps, num_of_frames * sizeof(uint64_t));
	memcpy(trajectory->initial_state, &initial_state, sizeof(NeedleCoordinates));
	memcpy(trajectory->final_state, &final_state, sizeof(NeedleCoordinates));

	if (PRINT_DEBUG_INFO) {
		print_trajectory(trajectory, num_of_frames, M_frames[num_of_frames-1]);
	}

	return trajectory;
}

// Python - C interface
// bool set_frames_py(
// 		double* table_x,
// 		double* table_y,
// 		double* table_z,
// 		double* slider_x,
// 		double* slider_y,
// 		double* slider_z,
// 		const double* ptr_state_A,
// 		const double* ptr_state_B,
// 		const size_t num_of_frames
// ){
// 	const NeedleCoordinates state_A = {ptr_state_A[0], ptr_state_A[1], ptr_state_A[2], ptr_state_A[3], ptr_state_A[4], ptr_state_A[5]};
// 	const NeedleCoordinates state_B = {ptr_state_B[0], ptr_state_B[1], ptr_state_B[2], ptr_state_B[3], ptr_state_B[4], ptr_state_B[5]};
//
// 	const Frames* frames = get_frames(state_A, state_B, num_of_frames);
//
// 	bool success = true;
// 	for (size_t i = 0; i < num_of_frames; ++i) {
// 		for (int j = 0; j < NUM_OF_VECTORS; ++j) {
// 			table_x[i * NUM_OF_VECTORS + j] = gsl_vector_get(frames->table[i].coordinates[j], 0);
// 			table_y[i * NUM_OF_VECTORS + j] = gsl_vector_get(frames->table[i].coordinates[j], 1);
// 			table_z[i * NUM_OF_VECTORS + j] = gsl_vector_get(frames->table[i].coordinates[j], 2);
// 		}
// 		for (int k = 0; k < NUM_OF_SLIDERS; ++k) {
// 			if (frames->sliders[i].coordinates[k] != NULL) {
// 				slider_x[i * NUM_OF_SLIDERS + k] = gsl_vector_get(frames->sliders[i].coordinates[k], 0);
// 				slider_y[i * NUM_OF_SLIDERS + k] = gsl_vector_get(frames->sliders[i].coordinates[k], 1);
// 				slider_z[i * NUM_OF_SLIDERS + k] = gsl_vector_get(frames->sliders[i].coordinates[k], 2);
// 			} else {
// 				slider_x[i * NUM_OF_SLIDERS + k] = NAN;
// 				slider_y[i * NUM_OF_SLIDERS + k] = NAN;
// 				slider_z[i * NUM_OF_SLIDERS + k] = NAN;
// 				success = false;
// 			}
// 		}
// 	}
// 	return success;
// }

uint64_t get_current_time_in_mks() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000 * 1000 + (uint64_t)ts.tv_nsec / 1000;
}

int compare_uint64_t(const void *a, const void *b) {
    uint64_t num1 = *(const uint64_t *)a;
    uint64_t num2 = *(const uint64_t *)b;

    if (num1 < num2) return -1;
    if (num1 > num2) return 1;
    return 0;
}


int main(void) {
	const NeedleCoordinates state_A = {-2.58, 1.77, 3.94, 0.71, 1.22, 0.47};
	const NeedleCoordinates state_B = {0.0, -1.45, 1.0, 0.20, -0.61, -0.41};
	const size_t num_of_frames = 10;

	uint64_t current_time = get_current_time_in_mks();
	Trajectory* trajectory = get_trajectory(state_A, state_B, current_time);
	
	uint64_t arbitrary_timestamp = current_time + (uint64_t)(1.2 * 1000 * 1000);
	printf("\nArbitrary timestamp:  %" PRIu64 "\n", arbitrary_timestamp);
	uint64_t *result = (uint64_t *)bsearch(&arbitrary_timestamp, trajectory->timestamps, num_of_frames, sizeof(uint64_t), compare_uint64_t);
	if (result != NULL) {
		printf("Exact timestamp found: %llu\n", *result);
	} else {
		int low = 0, high = num_of_frames - 1;

		while (low <= high) {
			int mid = (low + high) / 2;
			if (trajectory->timestamps[mid] < arbitrary_timestamp) {
				low = mid + 1;
			} else {
				high = mid - 1;
			}
		}

		uint64_t next_higher = (low < num_of_frames) ? trajectory->timestamps[low] : 0;  
		printf("Next timestamp found: %llu\n", next_higher);
	}
	return 0;
}
