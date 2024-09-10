#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // for memcpy()
#include <stdbool.h> // for bool
#include <time.h>   // for clock_gettime()
#include <unistd.h> // for usleep()
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#undef NDEBUG // Ensures that assertions are not disabled
#include <assert.h>


/// CONSTANTS THAT SHOULD NOT BE CHANGED
const double DEGREES_TO_RADIANS = (2 * M_PI / 360);
const size_t NUM_OF_HEXAGON_VERTICES = 6;
const size_t NUM_OF_VECTORS = NUM_OF_HEXAGON_VERTICES + 2;
const size_t NUM_OF_SLIDERS = 6;
///


/// CONSTANTS THAT CAN BE CHANGED
const bool PRINT_DEBUG_INFO_GLOBAL = true;
const size_t PRINT_DEBUG_SUBSELETION_FACTOR = 10 * 1000; // Prints info every nth iteration for some large logs;
		
const size_t NUM_OF_FRAMES = 10 * 1000;
const size_t MOVEMENT_DURATION_MKS = 5 * 1000 * 1000;
const uint64_t TRAJECTORY_WAIT_TIME_MKS = 500;
const double MIN_SLIDER_STEP_CM = 0.001;
const double NEEDLE_LENGTH_CM = 3.0;
const double TABLE_RADIUS_CM = 2.5; // Radius of the hexagonal table in cm.
const double TABLE_ANGLE = DEGREES_TO_RADIANS * 80; // Angle between vertices of the hexagonal table in radians. A regular hexagon will have 60 deg.
const double BASE_RADIUS_CM = 7.0; // Radius of the hexagonal base in cm.
const double BASE_ANGLE = DEGREES_TO_RADIANS * 15; // Angle between vertices of the hexagonal base in radians. A regular hexagon will have 60 deg.
const double MAX_POSSIBLE_HEIGHT_CM = 20.0; // Highest possible position for sliders in cm.
const double MIN_POSSIBLE_HEIGHT_CM = 0; // Lowest possible position for sliders in cm.
const double ARM_LENGTH_CM = 8.0; // Length of all robot arms in cm.
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

typedef struct {
	bool slider_movement_needed [NUM_OF_SLIDERS];
	bool slider_movement_direction_is_up [NUM_OF_SLIDERS];
} Sliders_Movement_Recipe;

typedef struct{
	// Pointers to arrays with size equal to the number of frames  
	Table_Model* table_frames; 
	Sliders_Model_Nullable* slider_frames;
	uint64_t* timestamps;

	// Pointers to single structs
	NeedleCoordinates* initial_state;
	NeedleCoordinates* final_state;

	// Number of frames
	size_t num_of_frames;
} Trajectory;

bool all_not_NULL(const Sliders_Model_Nullable sliders) {
    for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
        gsl_vector* vec_p = sliders.coordinates[i];
        if (vec_p == NULL) {
            return false;
        }
    }
    return true;
}
void print_needle_coordinates(NeedleCoordinates* state) {
	printf("Linear coordinates (x, y, z) in cm:        (%f, %f, %f)\n", state->x, state->y, state->z);
	printf("Euler angles (phi, theta, psi) in radians: (%f, %f, %f)\n", state->phi, state->theta, state->psi);
}

void print_table(const Table_Model table) {
	for(size_t i = 0; i < NUM_OF_VECTORS; i++) {

		printf("%lu: (%f, %f, %f)", i,
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
	for(size_t i = 0; i < NUM_OF_SLIDERS; i++) {
		if (sliders.coordinates[i] != NULL) {
			printf("Slider %lu: (%10.6f, %10.6f, %10.6f)\n", i,
					gsl_vector_get(sliders.coordinates[i], 0),
					gsl_vector_get(sliders.coordinates[i], 1),
					gsl_vector_get(sliders.coordinates[i], 2));
		} else {
			printf("Slider %lu: ERROR\n", i);
		}
	}
}

Sliders_Model_Nullable find_slider_difference(const Sliders_Model_Nullable sliders_2, const Sliders_Model_Nullable sliders_1) {
    Sliders_Model_Nullable difference;

    for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
        if (sliders_1.coordinates[i] != NULL && sliders_2.coordinates[i] != NULL) {
            difference.coordinates[i] = gsl_vector_alloc(3);
            
            for (size_t j = 0; j < 3; j++) {
                double diff = gsl_vector_get(sliders_2.coordinates[i], j) - gsl_vector_get(sliders_1.coordinates[i], j);
                gsl_vector_set(difference.coordinates[i], j, diff);
            }
        } else {
            difference.coordinates[i] = NULL;
        }
    }

    return difference;
}

void print_matrix_4D(const gsl_matrix* m) {
	for (size_t i = 0; i < 4; ++i)
		for (size_t j = 0; j < 4; ++j)
			printf(j==3?"%9.6f\n":"%9.6f ", gsl_matrix_get(m,i,j));
}

void print_trajectory(const Trajectory* trajectory) {
	const size_t num_of_frames = trajectory->num_of_frames;
	printf("\n==== Trajectory ====");
	printf("\nNumber of frames: %lu", num_of_frames);
	printf("\nTrajectory start time: %llu microseconds\n", trajectory->timestamps[0]);
	printf("Trajectory end time: %llu microseconds\n", trajectory->timestamps[num_of_frames - 1]);


	printf("\n==== Needle coordinates for state_initial:\n");
	print_needle_coordinates(trajectory->initial_state);
	printf("\n==== Needle coordinates for state_final:\n");
	print_needle_coordinates(trajectory->final_state);

	for (size_t i = 0; i < num_of_frames; i++) {
		bool print_debug_info_local = (i % PRINT_DEBUG_SUBSELETION_FACTOR) == 0;
		if (print_debug_info_local) {
			double progress = (double) i / (double) (num_of_frames - 1);
			printf("\n==== Frame %lu / %lu\nProgress: %f\nTimestamp: %llu \n", i, num_of_frames - 1, progress, trajectory->timestamps[i]);
			printf("\nTable coordinates (x, y, z):\n");
			print_table(trajectory->table_frames[i]);
			printf("\nSlider coordinates (x, y, z):\n");
			print_sliders(trajectory->slider_frames[i]);
		}
	}
}


void print_sliders_movement_recipe(const Sliders_Movement_Recipe* recipe) {
    printf("\nSliders Movement Recipe:\n");
	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
		printf("Slider %lu movement needed: %s\n", i, recipe->slider_movement_needed[i] ? "true" : "false");

	}
	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
		printf("Slider %lu movement direction is upwards: %s\n", i, recipe->slider_movement_direction_is_up[i] ? "true" : "false");
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
Sliders_Model_Nullable* get_sliders_in_state_0(const double radius, const double angle_in_radians, bool print_debug_info_local) {
	Sliders_Model_Nullable* sliders = (Sliders_Model_Nullable*) malloc(sizeof(Sliders_Model_Nullable));

	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
		sliders->coordinates[i] = gsl_vector_alloc(3);
	}

	gsl_vector_set(sliders->coordinates[0], 0, radius * sin(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[0], 1, radius * cos(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[0], 2, MAX_POSSIBLE_HEIGHT_CM);
	gsl_vector_set(sliders->coordinates[1], 0, radius * sin(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[1], 1, radius * cos(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[1], 2, MAX_POSSIBLE_HEIGHT_CM);
	gsl_vector_set(sliders->coordinates[2], 0, radius * sin(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[2], 1, radius * cos(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
gsl_vector_set(sliders->coordinates[2], 2, MAX_POSSIBLE_HEIGHT_CM);
	gsl_vector_set(sliders->coordinates[3], 0, radius * sin(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[3], 1, radius * cos(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[3], 2, MAX_POSSIBLE_HEIGHT_CM);
	gsl_vector_set(sliders->coordinates[4], 0, radius * sin(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[4], 1, radius * cos(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[4], 2, MAX_POSSIBLE_HEIGHT_CM);
	gsl_vector_set(sliders->coordinates[5], 0, radius * sin(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[5], 1, radius * cos(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[5], 2, MAX_POSSIBLE_HEIGHT_CM);


	if (print_debug_info_local) {
		printf("\n==== Slider coordinates (x, y, z) in cm for state_0:\n");
		print_sliders(*sliders);
	}

	return sliders;
}

Table_Model* get_table_in_state_0(bool print_debug_info_local) {
	Table_Model* table = (Table_Model*)malloc(sizeof(Table_Model));
	for (size_t i = 0; i < NUM_OF_VECTORS; i++) {
		table->coordinates[i] = gsl_vector_alloc(3);
	}

	Sliders_Model_Nullable* sliders_in_state_0 = get_sliders_in_state_0(BASE_RADIUS_CM, BASE_ANGLE, false);

	double radius = TABLE_RADIUS_CM;
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
	double deltaZ_squared = pow(ARM_LENGTH_CM, 2) - pow(deltaX, 2) - pow(deltaY, 2);
	assert(deltaZ_squared > 0);
	double table_z = slider_z - sqrt(deltaZ_squared); 

	for (size_t i = 0; i < NUM_OF_HEXAGON_VERTICES; i++) {
		gsl_vector_set(table->coordinates[i], 2, table_z);
	}

	gsl_vector_set(table->coordinates[6], 0, 0.0); 
	gsl_vector_set(table->coordinates[6], 1, 0.0);
	gsl_vector_set(table->coordinates[6], 2, table_z); 

	gsl_vector_set(table->coordinates[7], 0, 0.0);
	gsl_vector_set(table->coordinates[7], 1, 0.0);
	gsl_vector_set(table->coordinates[7], 2, table_z - NEEDLE_LENGTH_CM);

	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
		gsl_vector_free(sliders_in_state_0->coordinates[i]);
	}
	free(sliders_in_state_0);

	if (print_debug_info_local) {
		printf("\n==== Table coordinates (x, y, z) in cm for state_0:\n");
		print_table(*table);
	}

	return table;
}

Sliders_Model_Nullable* find_sliders(const Table_Model table_transformed_vecs, const Sliders_Model_Nullable sliders_in_state_0) {
	Sliders_Model_Nullable* sliders = (Sliders_Model_Nullable*) malloc(sizeof(Sliders_Model_Nullable));

	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
		double slider_x = gsl_vector_get(sliders_in_state_0.coordinates[i], 0);
		double slider_y = gsl_vector_get(sliders_in_state_0.coordinates[i], 1);

		double deltaX = gsl_vector_get(sliders_in_state_0.coordinates[i], 0) - gsl_vector_get(table_transformed_vecs.coordinates[i], 0);
		double deltaY = gsl_vector_get(sliders_in_state_0.coordinates[i], 1) - gsl_vector_get(table_transformed_vecs.coordinates[i], 1);
		double deltaZ_squared = pow(ARM_LENGTH_CM, 2) - pow(deltaX, 2) - pow(deltaY, 2);

		if (deltaZ_squared > 0) {
			// This is one of two possible solutions for z coordinate of slider, the higher one;
			// Lower solution will have minus sqrt(targerted_deltaZ_squared) instead of plus sqrt(targerted_deltaZ_squared);
			// We don't use second solution at alle because we fear damage during possible solution switching.
			double slider_z = gsl_vector_get(table_transformed_vecs.coordinates[i], 2) + sqrt(deltaZ_squared);
			if ((slider_z < MAX_POSSIBLE_HEIGHT_CM) && (slider_z > MIN_POSSIBLE_HEIGHT_CM)) {
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

	for (size_t i = 0; i < NUM_OF_VECTORS; i++) {
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

NeedleCoordinates* get_state_0(bool print_debug_info_local) {
	const Table_Model table_in_state_0 = *get_table_in_state_0(false);
	double needle_z_in_state_0 = gsl_vector_get(table_in_state_0.coordinates[7], 2);

	NeedleCoordinates* state_0 = (NeedleCoordinates*) malloc(sizeof(NeedleCoordinates));
	*state_0 = (NeedleCoordinates) {0.0, 0.0, needle_z_in_state_0, 0.0, 0.0, 0.0};

	if (print_debug_info_local) {
		printf("==== Needle coordinates for state_0:\n");
		print_needle_coordinates(state_0);
	}
	return state_0;
}

Trajectory* get_trajectory(const NeedleCoordinates state_initial, const NeedleCoordinates state_final, uint64_t movement_start_mks) {
	const uint64_t movement_duration_mks = MOVEMENT_DURATION_MKS; // TODO - calculate this
	const size_t num_of_frames = NUM_OF_FRAMES; // TODO - calculate this
	assert(num_of_frames >= 2);

	const NeedleCoordinates state_0 = *get_state_0(PRINT_DEBUG_INFO_GLOBAL);
	// Struct with array of pointers to vectors that describe Table 3D model in state_0.
	const Table_Model table_in_state_0 = *get_table_in_state_0(PRINT_DEBUG_INFO_GLOBAL);

	const gsl_matrix* M_0A = get_transformation_matrix(state_0, state_initial);

	// Struct with array of pointers to vectors that describe Table 3D model in state_A
	Table_Model table_in_state_A = *get_transformed_table(M_0A, table_in_state_0);

	// Struct with array of pointers to vectors that describe Sliders 3D model in zero state
	Sliders_Model_Nullable sliders_in_state_0 = *get_sliders_in_state_0(BASE_RADIUS_CM, BASE_ANGLE, PRINT_DEBUG_INFO_GLOBAL);

	// First frame is the the identity transformation (Id).
	// Last frame is the full transformation (M_AB), which brings (state_A) to the (state_B).
	// Other frames describe transformation that is "in progress", which bring (state_A) to some intermediate states (state_Tn)
	// Example for the case with 4 frames:
	// 	frame 0 is the Id;
	//	frame 1 is the M_AT1
	//	frame 2 is the M_AT2
	//	frame 3 is the M_AB;
	gsl_matrix* M_frames[num_of_frames];
	set_transformation_frames(M_frames, state_initial, state_final, num_of_frames);

	// Each 3D model is produced by applying one of transformations "in progress" to the 3D model for initial state.
	Table_Model table_frames[num_of_frames];
	set_table_frames(table_frames, num_of_frames, (const gsl_matrix**) M_frames, table_in_state_A);

	Sliders_Model_Nullable slider_frames[num_of_frames];
	set_sliders_frames(slider_frames, num_of_frames, table_frames, sliders_in_state_0);
	// print_transformed_slider_coordinates(slider_frames, num_of_frames);
	//
	uint64_t timestamps[num_of_frames];
	for (size_t i = 0; i < num_of_frames; i++) {
    		const double progress = (double) i / (double)(num_of_frames - 1);
	    	timestamps[i] = movement_start_mks + (uint64_t)(movement_duration_mks * progress);
	}

	Trajectory* trajectory = (Trajectory*)malloc(sizeof(Trajectory));
	trajectory->table_frames = (Table_Model*)malloc(num_of_frames * sizeof(Table_Model));
	trajectory->slider_frames = (Sliders_Model_Nullable*)malloc(num_of_frames * sizeof(Sliders_Model_Nullable));
	trajectory->timestamps = (uint64_t*)malloc(num_of_frames * sizeof(uint64_t));
	trajectory->initial_state = malloc(sizeof(NeedleCoordinates));
	trajectory->final_state = malloc(sizeof(NeedleCoordinates));

	memcpy(trajectory->table_frames, table_frames, num_of_frames * sizeof(Table_Model));
	memcpy(trajectory->slider_frames, slider_frames, num_of_frames * sizeof(Sliders_Model_Nullable));
	memcpy(trajectory->timestamps, timestamps, num_of_frames * sizeof(uint64_t));
	memcpy(trajectory->initial_state, &state_initial, sizeof(NeedleCoordinates));
	memcpy(trajectory->final_state, &state_final, sizeof(NeedleCoordinates));
	trajectory->num_of_frames = num_of_frames;

	if (PRINT_DEBUG_INFO_GLOBAL) {
		print_trajectory(trajectory);
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

size_t get_next_timestamp_index(uint64_t input_timestamp, size_t num_of_timestamps, uint64_t* timestamps, bool print_debug_info_local) {
    size_t result_index = SIZE_MAX; // This value indicates failure

    if (input_timestamp < timestamps[0]) {
	result_index = 0;
    }
    else if (input_timestamp > timestamps[num_of_timestamps - 1]) {
        result_index = SIZE_MAX;
    } else {
        // Binary search loop to find the exact or next bigger timestamp
        int moving_start_index = 0;
        int moving_end_index = num_of_timestamps - 1;

        while (moving_start_index <= moving_end_index) {
            int current_index = (moving_start_index + moving_end_index) / 2;
            if (timestamps[current_index] == input_timestamp) {
		result_index = current_index;
                break;
            } else if (timestamps[current_index] < input_timestamp) {
                moving_start_index = current_index + 1;
            } else {
                moving_end_index = current_index - 1;
            }
        }

        if (result_index == SIZE_MAX && moving_start_index < num_of_timestamps) {
  	    result_index = moving_start_index;
        }
    }

    if (PRINT_DEBUG_INFO_GLOBAL && print_debug_info_local) {
        if (result_index != SIZE_MAX) {
            printf("\nSearching for timestamp %llu\nNext timestamp found:   %llu (frame %lu / %lu)\n", input_timestamp, timestamps[result_index], result_index, num_of_timestamps - 1);
        } else {
            printf("\nSearching for timestamp %llu\nNo next timestamp found\n", input_timestamp);
        }
    }

    return result_index;
}

Sliders_Movement_Recipe compare_with_trajectory(const Trajectory* trajectory, const Sliders_Model_Nullable* sliders_inferred, const size_t frame_index, bool print_debug_info_local){
	Sliders_Model_Nullable sliders_at_trajectory = trajectory->slider_frames[frame_index];
	Sliders_Model_Nullable difference = find_slider_difference(sliders_at_trajectory, *sliders_inferred);
	Sliders_Movement_Recipe sliders_movement_recipe = {
		.slider_movement_needed[0] = fabs(gsl_vector_get(difference.coordinates[0], 2)) >= MIN_SLIDER_STEP_CM, 
		.slider_movement_needed[1] = fabs(gsl_vector_get(difference.coordinates[1], 2)) >= MIN_SLIDER_STEP_CM, 
		.slider_movement_needed[2] = fabs(gsl_vector_get(difference.coordinates[2], 2)) >= MIN_SLIDER_STEP_CM, 
		.slider_movement_needed[3] = fabs(gsl_vector_get(difference.coordinates[3], 2)) >= MIN_SLIDER_STEP_CM, 
		.slider_movement_needed[4] = fabs(gsl_vector_get(difference.coordinates[4], 2)) >= MIN_SLIDER_STEP_CM, 
		.slider_movement_needed[5] = fabs(gsl_vector_get(difference.coordinates[5], 2)) >= MIN_SLIDER_STEP_CM, 
		.slider_movement_direction_is_up[0] = gsl_vector_get(difference.coordinates[0], 2) >= 0.0,
		.slider_movement_direction_is_up[1] = gsl_vector_get(difference.coordinates[1], 2) >= 0.0,
		.slider_movement_direction_is_up[2] = gsl_vector_get(difference.coordinates[2], 2) >= 0.0,
		.slider_movement_direction_is_up[3] = gsl_vector_get(difference.coordinates[3], 2) >= 0.0,
		.slider_movement_direction_is_up[4] = gsl_vector_get(difference.coordinates[4], 2) >= 0.0,
		.slider_movement_direction_is_up[5] = gsl_vector_get(difference.coordinates[5], 2) >= 0.0,
	};

	if (PRINT_DEBUG_INFO_GLOBAL && print_debug_info_local) {
		printf("\n=== Comparing current sliders and sliders expected for frame (%lu / %lu)", frame_index, trajectory->num_of_frames -  1);
		printf("\nCurrent sliders (inferred from past inputs to the physical system):\n");
		print_sliders(*sliders_inferred);
		printf("\nExpected sliders:\n");
		print_sliders(sliders_at_trajectory);
		printf("\nDifference:\n");
		print_sliders(difference);
		print_sliders_movement_recipe(&sliders_movement_recipe);
	}

	return sliders_movement_recipe;
}

void mutate_infered_slider_state (Sliders_Model_Nullable* current_sliders, const Sliders_Movement_Recipe sliders_movement_recipe) {
	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
		if (sliders_movement_recipe.slider_movement_needed[i]) {
			if (sliders_movement_recipe.slider_movement_direction_is_up[i]) { 
				gsl_vector_set(current_sliders->coordinates[i], 2, gsl_vector_get(current_sliders->coordinates[i], 2) + MIN_SLIDER_STEP_CM);
			} else {
				gsl_vector_set(current_sliders->coordinates[i], 2, gsl_vector_get(current_sliders->coordinates[i], 2) - MIN_SLIDER_STEP_CM);
			}
		} 
	}
}

int follow_trajectory(Sliders_Model_Nullable* current_sliders, const Trajectory* trajectory) {
	const size_t num_of_frames = trajectory->num_of_frames;

	size_t iterations = 50 * 1000;
	size_t sleep_duration_mks = MOVEMENT_DURATION_MKS / iterations;

	for (size_t i = 0; i <= iterations; i++) {
		bool print_debug_info_local = (i % PRINT_DEBUG_SUBSELETION_FACTOR) == 0;
		size_t index_of_the_next_frame = get_next_timestamp_index(get_current_time_in_mks(), num_of_frames, trajectory->timestamps, print_debug_info_local);

		// if next frame found
		if (index_of_the_next_frame != SIZE_MAX) {
			Sliders_Movement_Recipe recipe = compare_with_trajectory(trajectory, current_sliders, index_of_the_next_frame, print_debug_info_local);
			mutate_infered_slider_state(current_sliders, recipe);
			if (
					!recipe.slider_movement_needed[0] &&
					!recipe.slider_movement_needed[1] &&
					!recipe.slider_movement_needed[2] &&
					!recipe.slider_movement_needed[3] &&
					!recipe.slider_movement_needed[4] &&
					index_of_the_next_frame == num_of_frames - 1
			   ) {
				if (PRINT_DEBUG_INFO_GLOBAL) {
					printf("\n=== Trajectory from (initial_state) to (final_state) followed sucessfully:\n");
					printf("\nstate_initial:\n");
					print_needle_coordinates(trajectory->initial_state);
					printf("\nstate_final:\n");
					print_needle_coordinates(trajectory->final_state);
				}
				return 0;
			}
		}
	usleep(sleep_duration_mks);
	} 

	return 1;
}

int main(void) {
	const NeedleCoordinates state_A = {-2.58, 1.77, 3.94, 0.71, 1.22, 0.47};
	const NeedleCoordinates state_B = {0.0, -1.45, 1.0, 0.20, -0.61, -0.41};

	const Trajectory* trajectory = get_trajectory(state_A, state_B, get_current_time_in_mks() + TRAJECTORY_WAIT_TIME_MKS);
	Sliders_Model_Nullable* current_sliders = malloc(sizeof(Sliders_Model_Nullable));
	current_sliders = &trajectory->slider_frames[0];
	
	(void) follow_trajectory(current_sliders, trajectory);

	return 0;
}
	
