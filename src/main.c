#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <string.h>

// Ensures that assertions are not disabled
#undef NDEBUG
#include <assert.h>


/// CONSTANTS THAT SHOULD NOT BE CHANGED
const double DEGREES_TO_RADIANS = (2 * M_PI / 360);
const int NUM_OF_HEXAGON_VERTICES = 6;
const int NUM_OF_VECTORS = NUM_OF_HEXAGON_VERTICES + 2;
const int NUM_OF_SLIDERS = 6;
///


/// CONSTANTS THAT CAN BE CHANGED
/*const int NUM_OF_HOMOTOPY_FRAMES = 3;*/
const double ACCURACY = 0.99; 
const double NEEDLE_LENGTH = 3;
const double TABLE_RADIUS = 3.0; // Radius of the hexagonal table in cm.
const double TABLE_ANGLE = DEGREES_TO_RADIANS * 80; // Angle between vertices of the hexagonal table in radians. A regular hexagon will have 60 deg.
const double BASE_RADIUS = 7.0; // Radius of the hexagonal base in cm.
const double BASE_ANGLE = DEGREES_TO_RADIANS * 15; // Angle between vertices of the hexagonal base in radians. A regular hexagon will have 60 deg.
const double MAX_POSSIBLE_HEIGHT = 15.0; // Highest possible position for sliders in cm.
const double MIN_POSSIBLE_HEIGHT = 0; // Lowest possible position for sliders in cm.
const double ARM_LENGTH = 8.0; // Length of all robot arms in cm.
///

typedef struct { // Full ideal description of the needle position, and thus, whole system.
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
	Table_Model* table; 
	Sliders_Model_Nullable* sliders; 
} Frames;

bool all_not_NULL(const Sliders_Model_Nullable sliders) {
    for (int i = 0; i < NUM_OF_SLIDERS; ++i) {
        gsl_vector* vec_p = sliders.coordinates[i];
        if (vec_p == NULL) {
            return false;
        }
    }
    return true;
}

void print_needle_coordinates(NeedleCoordinates table) {
    printf("Linear coordinates (x, y, z) in cm:        (%f, %f, %f)\n", table.x, table.y, table.z);
    printf("Euler angles (phi, theta, psi) in radians: (%f, %f, %f)\n", table.phi, table.theta, table.psi);
}

void print_matrix_4D(const gsl_matrix* m) {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            printf(j==3?"%9.6f\n":"%9.6f ", gsl_matrix_get(m,i,j));
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

void print_transformed_slider_coordinates(const Sliders_Model_Nullable* slider_frames, const size_t num_of_frames) {
	for (size_t j = 0; j < num_of_frames; j++)  {
		bool slider_coordinates_found = all_not_NULL(slider_frames[j]);
		if (slider_coordinates_found) {
			printf("\n==== Slider coordinates (x, y, z) in cm for frame %d:\n", (int) j );
			print_sliders(slider_frames[j]);

		} else {
			printf("\n==== Sliders coordinates for frame %d NOT found!\n", (int) j);
			printf("Failed to find slider coordinates that can realize targeted table coordinates.\n");
			printf("Targeted table coordinates deemed infeasible. Slider diagnostics:\n");
			print_sliders(slider_frames[j]);
		}
	}
}

gsl_matrix* get_transformation_matrix(
    const NeedleCoordinates state_A,
    const NeedleCoordinates state_B
    ) {
	printf("\n==== Transformation matrices for transition between state_A and state_B");
	printf("\n==== state_A:\n");
	print_needle_coordinates(state_A);
	printf("\n==== state_B:\n");
	print_needle_coordinates(state_B);

	gsl_matrix* T_0B = gsl_matrix_alloc(4, 4);
	gsl_matrix_set(T_0B, 0, 0, cos(state_B.psi) * cos(state_B.theta));
	gsl_matrix_set(T_0B, 0, 1, cos(state_B.psi) * sin(state_B.theta) * sin(state_B.phi) - sin(state_B.psi) * cos(state_B.phi));
	gsl_matrix_set(T_0B, 0, 2, cos(state_B.psi) * sin(state_B.theta) * cos(state_B.phi) + sin(state_B.psi) * sin(state_B.phi));
	gsl_matrix_set(T_0B, 0, 3, state_B.x);
	gsl_matrix_set(T_0B, 1, 0, sin(state_B.psi) * cos(state_B.theta));
	gsl_matrix_set(T_0B, 1, 1, sin(state_B.psi) * sin(state_B.theta) * sin(state_B.phi) + cos(state_B.psi) * cos(state_B.phi));
	gsl_matrix_set(T_0B, 1, 2, sin(state_B.psi) * sin(state_B.theta) * cos(state_B.phi) - cos(state_B.psi) * sin(state_B.phi));
	gsl_matrix_set(T_0B, 1, 3, state_B.y);
	gsl_matrix_set(T_0B, 2, 0, -sin(state_B.theta));
	gsl_matrix_set(T_0B, 2, 1, cos(state_B.theta) * sin(state_B.phi));
	gsl_matrix_set(T_0B, 2, 2, cos(state_B.theta) * cos(state_B.phi));
	gsl_matrix_set(T_0B, 2, 3, state_B.z);
	gsl_matrix_set(T_0B, 3, 0, 0.0);
	gsl_matrix_set(T_0B, 3, 1, 0.0);
	gsl_matrix_set(T_0B, 3, 2, 0.0);
	gsl_matrix_set(T_0B, 3, 3, 1.0);
	printf("\n==== T_0B, transformation matrix for (0,0,0,0,0,0) -> (state_B):\n");
	print_matrix_4D(T_0B);

	gsl_matrix* T_0A = gsl_matrix_alloc(4, 4);
	gsl_matrix_set(T_0A, 0, 0, cos(state_A.psi) * cos(state_A.theta));
	gsl_matrix_set(T_0A, 0, 1, cos(state_A.psi) * sin(state_A.theta) * sin(state_A.phi) - sin(state_A.psi) * cos(state_A.phi));
	gsl_matrix_set(T_0A, 0, 2, cos(state_A.psi) * sin(state_A.theta) * cos(state_A.phi) + sin(state_A.psi) * sin(state_A.phi));
	gsl_matrix_set(T_0A, 0, 3, state_A.x);
	gsl_matrix_set(T_0A, 0, 3, state_A.x);
	gsl_matrix_set(T_0A, 1, 0, sin(state_A.psi) * cos(state_A.theta));
	gsl_matrix_set(T_0A, 1, 1, sin(state_A.psi) * sin(state_A.theta) * sin(state_A.phi) + cos(state_A.psi) * cos(state_A.phi));
	gsl_matrix_set(T_0A, 1, 2, sin(state_A.psi) * sin(state_A.theta) * cos(state_A.phi) - cos(state_A.psi) * sin(state_A.phi));
	gsl_matrix_set(T_0A, 1, 3, state_A.y);
	gsl_matrix_set(T_0A, 2, 0, -sin(state_A.theta));
	gsl_matrix_set(T_0A, 2, 1, cos(state_A.theta) * sin(state_A.phi));
	gsl_matrix_set(T_0A, 2, 2, cos(state_A.theta) * cos(state_A.phi));
	gsl_matrix_set(T_0A, 2, 3, state_A.z);
	gsl_matrix_set(T_0A, 3, 0, 0.0);
	gsl_matrix_set(T_0A, 3, 1, 0.0);
	gsl_matrix_set(T_0A, 3, 2, 0.0);
	gsl_matrix_set(T_0A, 3, 3, 1.0);
	printf("\n==== T_0A, transformation matrix for (0,0,0,0,0,0) -> (state_A):\n");
	print_matrix_4D(T_0A);

	int signum;
	gsl_permutation * p = gsl_permutation_alloc(4);
	gsl_matrix* T_A0 = gsl_matrix_alloc(4, 4);
	gsl_matrix* T_0A_temp = gsl_matrix_alloc(4, 4);
	gsl_matrix_memcpy(T_0A_temp, T_0A);
	gsl_linalg_LU_decomp (T_0A_temp, p, &signum);
	gsl_linalg_LU_invert (T_0A_temp, p, T_A0);
	printf("\n==== T_A0, transformation matrix for (state_A) -> (0,0,0,0,0,0), an inverse of T_0A:\n");
	print_matrix_4D(T_A0);
	gsl_permutation_free(p);

	// T_AB = T_0B * (T_0A)^(-1)
	// T_AB = T_0B * T_A0
	gsl_matrix* T_AB = gsl_matrix_alloc(4, 4);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T_0B, T_A0,  0.0, T_AB);
	printf("\n==== T_AB, transformation matrix for (state_A) -> (state_B):\n");
	print_matrix_4D(T_AB);

	return T_AB;
}

void set_transformation_frames(gsl_matrix** transformation_frames, const gsl_matrix* T, const size_t num_of_frames) {
	printf("\n==== Homotopy frames"); 

	for (size_t j = 0; j < num_of_frames; j++) {
		transformation_frames[j] = gsl_matrix_alloc(4, 4);
		gsl_matrix_memcpy(transformation_frames[j], T);

		// Here we mutate affine transformation matrix T into
		// (k)T + (1-k)Id, where k is the homotopy coefficient
		double homotopy_coefficient = (double) j / (double) (num_of_frames-1);
		gsl_matrix_scale(transformation_frames[j], homotopy_coefficient);
		gsl_matrix* Id = gsl_matrix_alloc(4, 4);
		gsl_matrix_set_identity(Id);
		gsl_matrix_scale(Id, 1.0 - homotopy_coefficient);
		gsl_matrix_add(transformation_frames[j], Id);

		printf("\n==== Transformation matrix for frame %d: (%.3f)T + (%.3f)Id:\n", (int) j, (homotopy_coefficient), 1-homotopy_coefficient);
		print_matrix_4D(transformation_frames[j]);
	}
}


Table_Model* get_model_in_state_0 () {
	Table_Model* model = (Table_Model*)malloc(sizeof(Table_Model));
	for (int i = 0; i < NUM_OF_VECTORS; i++) {
		model->coordinates[i] = gsl_vector_alloc(3);
	}

	double radius = TABLE_RADIUS;
	double angle_in_radians = TABLE_ANGLE;
	double distance_from_needle_end = NEEDLE_LENGTH;

	gsl_vector_set(model->coordinates[0], 0, radius * sin(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(model->coordinates[0], 1, radius * cos(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(model->coordinates[0], 2, distance_from_needle_end);
	gsl_vector_set(model->coordinates[1], 0, radius * sin(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(model->coordinates[1], 1, radius * cos(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(model->coordinates[1], 2, distance_from_needle_end);
	gsl_vector_set(model->coordinates[2], 0, radius * sin(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(model->coordinates[2], 1, radius * cos(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(model->coordinates[2], 2, distance_from_needle_end);
	gsl_vector_set(model->coordinates[3], 0, radius * sin(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(model->coordinates[3], 1, radius * cos(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(model->coordinates[3], 2, distance_from_needle_end);
	gsl_vector_set(model->coordinates[4], 0, radius * sin(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(model->coordinates[4], 1, radius * cos(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(model->coordinates[4], 2, distance_from_needle_end);
	gsl_vector_set(model->coordinates[5], 0, radius * sin(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(model->coordinates[5], 1, radius * cos(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(model->coordinates[5], 2, distance_from_needle_end);

	// Table centroid
	gsl_vector_set(model->coordinates[6], 0, 0.0);
	gsl_vector_set(model->coordinates[6], 1, 0.0);
	gsl_vector_set(model->coordinates[6], 2, distance_from_needle_end);

	// Needle end
	gsl_vector_set(model->coordinates[7], 0, 0.0);
	gsl_vector_set(model->coordinates[7], 1, 0.0);
	gsl_vector_set(model->coordinates[7], 2, 0.0);

	printf("\n==== Table point coordinates (x, y, z) in cm for state_0:\n");
	print_table(*model);

	return model;
}

// Coordinates of sliders at their lowest possible position at rest.
Sliders_Model_Nullable* get_sliders_in_state_0( const double radius, const double angle_in_radians) {
	Sliders_Model_Nullable* sliders = (Sliders_Model_Nullable*) malloc(sizeof(Sliders_Model_Nullable));

	for (int i = 0; i < NUM_OF_SLIDERS; i++) {
		sliders->coordinates[i] = gsl_vector_alloc(3);
	}

	gsl_vector_set(sliders->coordinates[0], 0, radius * sin(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[0], 1, radius * cos(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[0], 2, 0.0);
	gsl_vector_set(sliders->coordinates[1], 0, radius * sin(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[1], 1, radius * cos(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[1], 2, 0.0);
	gsl_vector_set(sliders->coordinates[2], 0, radius * sin(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[2], 1, radius * cos(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[2], 2, 0.0);
	gsl_vector_set(sliders->coordinates[3], 0, radius * sin(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[3], 1, radius * cos(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[3], 2, 0.0);
	gsl_vector_set(sliders->coordinates[4], 0, radius * sin(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[4], 1, radius * cos(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[4], 2, 0.0);
	gsl_vector_set(sliders->coordinates[5], 0, radius * sin(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[5], 1, radius * cos(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->coordinates[5], 2, 0.0);

	printf("\n==== Slider coordinates (x, y, z) in cm for state_0:\n");
	
	print_sliders(*sliders);
	return sliders;
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



double get_velocity(const double progress, const double velocity_max) {
	// Logistic curve is used:
	// https://www.desmos.com/calculator/j0sgzmoqbj
	assert(velocity_max > 0.0);
	assert(ACCURACY > 0.0 && ACCURACY < 1.0);
	assert(progress >= 0.0 && progress <= 1.0);

	double k = log((1.0-ACCURACY)/(ACCURACY));
	double velocity = velocity_max/(1.0 + exp(k*(2.0*progress - 1.0)));
	return velocity;
}

double get_acceleration(const double progress, const double velocity_max) {
	// Logistic curve is used:
	// https://www.desmos.com/calculator/j0sgzmoqbj
	assert(velocity_max > 0.0);
	assert(ACCURACY > 0.0 && ACCURACY < 1.0);
	assert(progress >= 0.0 && progress <= 1.0);

	double k = log((1.0-ACCURACY)/(ACCURACY));
	double acceleration = -2*k*velocity_max * exp(k*(2.0*progress - 1.0)) / pow(exp(k*(2.0*progress - 1.0))+1, 2);
	return acceleration;
}

Table_Model* get_transformed_model(const gsl_matrix* T, const Table_Model initial_model){
	Table_Model* transformed_model = (Table_Model*)malloc(sizeof(Table_Model));

	for (int i = 0; i < NUM_OF_VECTORS; i++) {
		const gsl_vector* initial_vec_in_3D = initial_model.coordinates[i];

		// Uplifting to 4D
		gsl_vector* initial_vec_in_4D = gsl_vector_alloc(4);
		gsl_vector_set(initial_vec_in_4D, 0, gsl_vector_get(initial_vec_in_3D, 0)); // We take X coordinate of 3D vector
		gsl_vector_set(initial_vec_in_4D, 1, gsl_vector_get(initial_vec_in_3D, 1)); // We take Y coordinate of 3D vector
		gsl_vector_set(initial_vec_in_4D, 2, gsl_vector_get(initial_vec_in_3D, 2)); // We take Z coordinate of 3D vector
		gsl_vector_set(initial_vec_in_4D, 3, 1.0); // We use new coordinate for Tau

		// Linear transformation in 4D
		gsl_vector* transformed_vec_in_4D = gsl_vector_alloc(4);
		gsl_blas_dgemv(CblasNoTrans, 1.0, T, initial_vec_in_4D, 0.0, transformed_vec_in_4D);

		/*// Projection back to 3D*/
		transformed_model->coordinates[i] = gsl_vector_alloc(3);
		gsl_vector_set(transformed_model->coordinates[i], 0, gsl_vector_get(transformed_vec_in_4D, 0));
		gsl_vector_set(transformed_model->coordinates[i], 1, gsl_vector_get(transformed_vec_in_4D, 1));
		gsl_vector_set(transformed_model->coordinates[i], 2, gsl_vector_get(transformed_vec_in_4D, 2));

		gsl_vector_free(initial_vec_in_4D);
		gsl_vector_free(transformed_vec_in_4D);
	}
	return transformed_model;
}

void set_model_frames(Table_Model* model_frames, const size_t num_of_frames, const gsl_matrix** T_frames, const Table_Model initial_model) {
	for (size_t j = 0; j < num_of_frames; j++) {
		model_frames[j] = *get_transformed_model(T_frames[j], initial_model);
		printf("\n==== Table point coordinates (x, y, z) in cm for frame %d:\n", (int) j);
		print_table(model_frames[j]);
		// const gsl_matrix** prohibits this:
		// gsl_matrix_set(T_frames[j], 1, 1, 0.0);
	}
}

/// C-Python interface


Frames* get_frames(NeedleCoordinates state_A, NeedleCoordinates state_B, const size_t num_of_frames) {
	assert(num_of_frames >= 2);
	const NeedleCoordinates state_0 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	const gsl_matrix* T_0A = get_transformation_matrix(state_0, state_A);
	const gsl_matrix* T_AB = get_transformation_matrix(state_A, state_B);
	
	// Struct with array of pointers to vectors that describe Table 3D model in state_0 ({0, 0, 0, 0, 0, 0} state)
	Table_Model model_in_state_0 = *get_model_in_state_0();

	// Struct with array of pointers to vectors that describe Table 3D model in state_A
	Table_Model model_in_state_A = *get_transformed_model(T_0A, model_in_state_0);

	// Struct with array of pointers to vectors that describe Sliders 3D model in zero state ({0, 0, 0, 0, 0, 0} state)
	Sliders_Model_Nullable sliders_in_state_0 = *get_sliders_in_state_0(BASE_RADIUS, BASE_ANGLE);

	// Transformation frames are created using homotopy Id->T
	// First frame is the the identity transformation (Id).
	// Last frame is the full transformation (T_AB), which brings (state_A) to the (state_B).
	// Other frames are created by mixing T_AB and Id, those frames describe transformation that is "in progress".
	// Example for the case with 3 frames:
	// 	frame 0 is the Id;
	//	frame 1 is 0.5 Id + 0.5 T_AB;
	//	frame 2 is the T_AB;
	gsl_matrix* T_frames[num_of_frames];
	set_transformation_frames(T_frames, T_AB, num_of_frames);

	// Each 3D model is produced by applying one of transformations "in progress" to the 3D model for initial state.
	Table_Model table_frames[num_of_frames];
	set_model_frames(table_frames, num_of_frames, (const gsl_matrix**) T_frames, model_in_state_A);

	Sliders_Model_Nullable slider_frames[num_of_frames];
	set_sliders_frames(slider_frames, num_of_frames, table_frames, sliders_in_state_0);
	print_transformed_slider_coordinates(slider_frames, num_of_frames);

	Frames* f = (Frames*)malloc(sizeof(Frames));
	f->table = (Table_Model*)malloc(num_of_frames * sizeof(Table_Model));
	f->sliders = (Sliders_Model_Nullable*)malloc(num_of_frames * sizeof(Sliders_Model_Nullable));
	memcpy(f->table, table_frames, num_of_frames * sizeof(Table_Model));
	memcpy(f->sliders, slider_frames, num_of_frames * sizeof(Sliders_Model_Nullable));

	return f;
}

// Python - C interface
bool set_frames_py(
		double* table_x,
		double* table_y,
		double* table_z,
		double* slider_x,
		double* slider_y,
		double* slider_z,
		const double* state_1,
		const double* state_2,
		const size_t num_of_frames
){
	const NeedleCoordinates state_A = {state_1[0], state_1[1], state_1[2], state_1[3], state_1[4], state_1[5]};
	const NeedleCoordinates state_B = {state_2[0], state_2[1], state_2[2], state_2[3], state_2[4], state_2[5]};

	const Frames* frames = get_frames(state_A, state_B, num_of_frames);
	
	bool success = true;
	for (size_t i = 0; i < num_of_frames; ++i) {
		for (int j = 0; j < NUM_OF_VECTORS; ++j) {
			table_x[i * NUM_OF_VECTORS + j] = gsl_vector_get(frames->table[i].coordinates[j], 0);
			table_y[i * NUM_OF_VECTORS + j] = gsl_vector_get(frames->table[i].coordinates[j], 1);
			table_z[i * NUM_OF_VECTORS + j] = gsl_vector_get(frames->table[i].coordinates[j], 2);
		}
		for (int k = 0; k < NUM_OF_SLIDERS; ++k) {
			if (frames->sliders[i].coordinates[k] != NULL) {
				slider_x[i * NUM_OF_SLIDERS + k] = gsl_vector_get(frames->sliders[i].coordinates[k], 0);
				slider_y[i * NUM_OF_SLIDERS + k] = gsl_vector_get(frames->sliders[i].coordinates[k], 1);
				slider_z[i * NUM_OF_SLIDERS + k] = gsl_vector_get(frames->sliders[i].coordinates[k], 2);
			} else {
				slider_x[i * NUM_OF_SLIDERS + k] = NAN;
				slider_y[i * NUM_OF_SLIDERS + k] = NAN;
				slider_z[i * NUM_OF_SLIDERS + k] = NAN;
				success = false;
			}
		}
	}
	return success;
}


int main(void) {
	const NeedleCoordinates state_A = {-2.58, 1.77, 1.94, 0.71, 1.22, 0.47};
	const NeedleCoordinates state_B = {0.0, -1.45, -2.58, 0.20, -0.61, -0.41};
	const size_t num_of_frames = 10;

	(void) get_frames(state_A, state_B, num_of_frames);

	return 0;
}
