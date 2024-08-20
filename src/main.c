#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

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


typedef gsl_vector* Model_3D[NUM_OF_VECTORS];

struct NeedleCoordinates { // Full ideal description of the needle position, and thus, whole system.
    double x, y, z; // Coordinates of the needle end in cm.
    double phi, theta, psi; // Euler angles for the needle in radians.
};

// Array of pointers to NULL (if slider coordinates were not found)
// or pointers to slider coordinates (if slider coordinates were found).
typedef struct {
    gsl_vector* coordinates[NUM_OF_SLIDERS];
} OptionSliderCoordinates;

bool all_not_NULL(const OptionSliderCoordinates* options) {
    for (int i = 0; i < NUM_OF_SLIDERS; ++i) {
        gsl_vector* vec_p = options->coordinates[i];
        if (vec_p == NULL) {
            return false;
        }
    }
    return true;
}

void print_needle_coordinates(struct NeedleCoordinates table) {
    printf("Linear coordinates (x, y, z) in cm:        (%f, %f, %f)\n", table.x, table.y, table.z);
    printf("Euler angles (phi, theta, psi) in radians: (%f, %f, %f)\n", table.phi, table.theta, table.psi);
}

void print_matrix_4D(const gsl_matrix* m) {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            printf(j==3?"%9.6f\n":"%9.6f ", gsl_matrix_get(m,i,j));
}

void print_model(gsl_vector* vectors[NUM_OF_VECTORS]) {
    for(int i = 0; i < NUM_OF_VECTORS; i++) {

         printf("%d: (%f, %f, %f)", i,
             gsl_vector_get(vectors[i], 0),
             gsl_vector_get(vectors[i], 1),
             gsl_vector_get(vectors[i], 2));

        if (i == NUM_OF_VECTORS - 2) {
            printf(" <-- table centroid");
        }
        if (i == NUM_OF_VECTORS - 1){
            printf(" <-- needle end");
        }
        printf("\n");
    }
}

void print_lowest_possible_slider_coordinates(gsl_vector** vectors, int num_of_vectors_to_print) {
    assert(num_of_vectors_to_print == NUM_OF_SLIDERS);

    for(int i = 0; i < num_of_vectors_to_print; i++) {
        printf("Slider %d: (%f, %f, %f)\n", i,
            gsl_vector_get(vectors[i], 0),
            gsl_vector_get(vectors[i], 1),
            gsl_vector_get(vectors[i], 2));
    }
}

void print_transformed_slider_coordinates(const OptionSliderCoordinates option) {
    bool slider_coordinates_found = all_not_NULL(&option);

    if (slider_coordinates_found) {
        printf("\n=== Solution found! ===\n");
        printf("Targeted table coordinates can be realized using these slider coordinates (x, y, z) in cm:\n");
        for (int i = 0; i < NUM_OF_SLIDERS; i++) {
            printf("Slider %d: (%f, %f, %f)\n", i,
                gsl_vector_get(option.coordinates[i], 0),
                gsl_vector_get(option.coordinates[i], 1),
                gsl_vector_get(option.coordinates[i], 2));
        }
    } else {
        printf("\n=== Solution NOT found! ===\n");
        printf("Failed to find slider coordinates that can realize targeted table coordinates.\n");
        printf("Targeted table coordinates deemed infeasible. Slider diagnostics:\n");
        for (int i = 0; i < NUM_OF_SLIDERS; i++) {
            bool slider_is_ok = (option.coordinates[i] != NULL);
            if (slider_is_ok) {
                printf("Slider %d: ok \n", i);
            } else {
                printf("Slider %d: ERROR\n", i);
            }
        }
    }
    printf("\n");
}

void set_transformation_matrix(
    gsl_matrix* T_AB, // Place where final affine transformation matrix will be written
    const struct NeedleCoordinates state_A,
    const struct NeedleCoordinates state_B
    ) {
	printf("\n────────── Transformation matrices ──────────");
	printf("\nstate_A:\n");
	print_needle_coordinates(state_A);
	printf("\nstate_B:\n");
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
	printf("\nT_0B, transformation matrix for (0,0,0,0,0,0) -> (state_B)\n");
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
	printf("\nT_0A, transformation matrix for (0,0,0,0,0,0) -> (state_A)\n");
	print_matrix_4D(T_0A);

	int signum;
	gsl_permutation * p = gsl_permutation_alloc(4);
	gsl_matrix* T_A0 = gsl_matrix_alloc(4, 4);
	gsl_matrix* T_0A_temp = gsl_matrix_alloc(4, 4);
	gsl_matrix_memcpy(T_0A_temp, T_0A);
	gsl_linalg_LU_decomp (T_0A_temp, p, &signum);
	gsl_linalg_LU_invert (T_0A_temp, p, T_A0);
	printf("\nT_A0, transformation matrix for (state_A) -> (0,0,0,0,0,0), an inverse of T_0A\n");
	print_matrix_4D(T_A0);
	gsl_permutation_free(p);

	// T_AB = T_0B * (T_0A)^(-1)
	// T_AB = T_0B * T_A0
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T_0B, T_A0,  0.0, T_AB);
	printf("\nT_AB, transformation matrix for (state_A) -> (state_B)\n");
	print_matrix_4D(T_AB);
}

void set_sequence(gsl_matrix** transformation_frames, const gsl_matrix* T, const size_t num_of_frames) {
	printf("\n────────── Homotopy frames ──────────");
	for (int i = 0; i < num_of_frames; i++) {
		double homotopy_coefficient = (double) i / (double) (num_of_frames-1);
		transformation_frames[i] = gsl_matrix_alloc(4, 4);
		gsl_matrix *T_copy = transformation_frames[i];
		gsl_matrix_memcpy(T_copy, T);

		// Here we mutate affine transformation matrix T into
		// (k)T + (1-k)Id, where k is the homotopy coefficient
		gsl_matrix_scale(T_copy, homotopy_coefficient);
		gsl_matrix* Id = gsl_matrix_alloc(4, 4);
		gsl_matrix_set_identity(Id);
		gsl_matrix_scale(Id, 1.0 - homotopy_coefficient);
		gsl_matrix_add(T_copy, Id);

		printf("\n(%.3f)T + (%.3f)Id\n", (homotopy_coefficient), 1-homotopy_coefficient);
		print_matrix_4D(T_copy);
	}
}


void set_model_in_state_0 (gsl_vector** vecs) {
	for (int i = 0; i < NUM_OF_VECTORS; i++) {
		vecs[i] = gsl_vector_alloc(3);
	}

	double radius = TABLE_RADIUS;
	double angle_in_radians = TABLE_ANGLE;
	double distance_from_needle_end = NEEDLE_LENGTH;

	gsl_vector_set(vecs[0], 0, radius * sin(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(vecs[0], 1, radius * cos(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(vecs[0], 2, distance_from_needle_end);

	gsl_vector_set(vecs[1], 0, radius * sin(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(vecs[1], 1, radius * cos(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(vecs[1], 2, distance_from_needle_end);

	gsl_vector_set(vecs[2], 0, radius * sin(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(vecs[2], 1, radius * cos(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(vecs[2], 2, distance_from_needle_end);

	gsl_vector_set(vecs[3], 0, radius * sin(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(vecs[3], 1, radius * cos(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(vecs[3], 2, distance_from_needle_end);

	gsl_vector_set(vecs[4], 0, radius * sin(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(vecs[4], 1, radius * cos(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(vecs[4], 2, distance_from_needle_end);

	gsl_vector_set(vecs[5], 0, radius * sin(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(vecs[5], 1, radius * cos(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(vecs[5], 2, distance_from_needle_end);

	// Table centroid
	gsl_vector_set(vecs[6], 0, 0.0);
	gsl_vector_set(vecs[6], 1, 0.0);
	gsl_vector_set(vecs[6], 2, distance_from_needle_end);

	// Needle end
	gsl_vector_set(vecs[7], 0, 0.0);
	gsl_vector_set(vecs[7], 1, 0.0);
	gsl_vector_set(vecs[7], 2, 0.0);

	printf("\nInitial (x, y, z) coordinates of 3d model points in cm:\n");
	print_model(vecs);
}

// Used to populate slider base coordinates.
// Base coordinates are coordinates of sliders at their lowest possible position at rest.
void set_hexagon(
    gsl_vector** vecs,
    const double radius,
    const double angle_in_radians) {
    gsl_vector_set(vecs[0], 0, radius * sin(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
    gsl_vector_set(vecs[0], 1, radius * cos(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
    gsl_vector_set(vecs[0], 2, 0.0);

    gsl_vector_set(vecs[1], 0, radius * sin(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
    gsl_vector_set(vecs[1], 1, radius * cos(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
    gsl_vector_set(vecs[1], 2, 0.0);

    gsl_vector_set(vecs[2], 0, radius * sin(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
    gsl_vector_set(vecs[2], 1, radius * cos(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
    gsl_vector_set(vecs[2], 2, 0.0);

    gsl_vector_set(vecs[3], 0, radius * sin(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
    gsl_vector_set(vecs[3], 1, radius * cos(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
    gsl_vector_set(vecs[3], 2, 0.0);

    gsl_vector_set(vecs[4], 0, radius * sin(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
    gsl_vector_set(vecs[4], 1, radius * cos(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
    gsl_vector_set(vecs[4], 2, 0.0);

    gsl_vector_set(vecs[5], 0, radius * sin(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
    gsl_vector_set(vecs[5], 1, radius * cos(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
    gsl_vector_set(vecs[5], 2, 0.0);
}



OptionSliderCoordinates find_slider_vecs(
    gsl_vector** table_transformed_vecs,
    gsl_vector** base_vecs) {

    OptionSliderCoordinates maybe_slider_vecs;

    for (int i = 0; i < NUM_OF_SLIDERS; i++) {
        maybe_slider_vecs.coordinates[i] = gsl_vector_alloc(3);
        gsl_vector_memcpy(maybe_slider_vecs.coordinates[i], base_vecs[i]);

        double table_transformed_x = gsl_vector_get(table_transformed_vecs[i], 0);
        double table_transformed_y = gsl_vector_get(table_transformed_vecs[i], 1);
        double table_transformed_z = gsl_vector_get(table_transformed_vecs[i], 2);
        double base_x = gsl_vector_get(base_vecs[i], 0);
        double base_y = gsl_vector_get(base_vecs[i], 1);

        bool slider_is_inside_range = false;
        // deltaZ^2 = R^2 - deltaX^2 - deltaY^2
        double z_squared = pow(ARM_LENGTH, 2) - pow((base_x - table_transformed_x), 2) - pow((base_y - table_transformed_y), 2);

        if (z_squared > 0) {
            // This is one of two possible solutions for z coordinate of slider (higher solution of two).
            // Lower solution will be table_transformed_z - sqrt(z_squared).
            // We don't use second solution because we fear damage during solution switching.
            double slider_z = table_transformed_z + sqrt(z_squared);
            slider_is_inside_range = (slider_z < MAX_POSSIBLE_HEIGHT) && (slider_z > MIN_POSSIBLE_HEIGHT);
            if (slider_is_inside_range) {
                gsl_vector_set(maybe_slider_vecs.coordinates[i], 2, slider_z);
            } else {
                gsl_vector_free(maybe_slider_vecs.coordinates[i]);
                maybe_slider_vecs.coordinates[i] = NULL;
            };
        } else {
            gsl_vector_free(maybe_slider_vecs.coordinates[i]);
            maybe_slider_vecs.coordinates[i] = NULL;
        }
    }
    return maybe_slider_vecs;
}

//typedef gsl_m?atrix gsl_matrix;

/// C-Python interface
/* double** set_transformed_model(const double* state_1, const double* state_2, const double* model_for_state_1, const size_t num_of_homotopy_frames){*/
/*     struct NeedleCoordinates state_1_struct = {state_1[0], state_1[1], state_1[2], state_1[3], state_1[4], state_1[5]};*/
/*     struct NeedleCoordinates state_2_struct = {state_2[0], state_2[1], state_2[2], state_2[3], state_2[4], state_2[5]};*/
/**/
/*     gsl_matrix* T_1to2 = gsl_matrix_alloc(4, 4);*/
/*     set_transformation_matrix(T_1to2, state_1_struct, state_2_struct);*/
/**/
/**/
/**/
/**/
/*     return model_for_state_2;*/
/* }*/
/*///*/


double targeted_velocity(const double progress, const double velocity_max) {
	// Logistic curve is used:
	// https://www.desmos.com/calculator/z8iodbso4f
	assert(velocity_max > 0.0);
	assert(ACCURACY > 0.0 && ACCURACY < 1.0);
	assert(progress >= 0.0 && progress <= 1.0);

	double k = log((1.0-ACCURACY)/(ACCURACY));
	double velocity = velocity_max/(1.0 + exp(k*(2.0*progress - 1.0)));
	return velocity;
}

double targeted_acceleration(const double progress, const double velocity_max) {
	// https://www.desmos.com/calculator/j0sgzmoqbj
	assert(velocity_max > 0.0);
	assert(ACCURACY > 0.0 && ACCURACY < 1.0);
	assert(progress >= 0.0 && progress <= 1.0);

	double k = log((1.0-ACCURACY)/(ACCURACY));
	double acceleration = -2*k*velocity_max * exp(k*(2.0*progress - 1.0)) / pow(exp(k*(2.0*progress - 1.0))+1, 2);
	return acceleration;
}

void set_transformed_model(Model_3D place_for_transformed_model, const gsl_matrix* T, const Model_3D initial_model ){
	for (int i = 0; i < NUM_OF_VECTORS; i++) {
		place_for_transformed_model[i] = gsl_vector_alloc(3);
		const gsl_vector* initial_vec_in_3D = initial_model[i];
		gsl_vector* place_for_transformed_vec_in_3D = place_for_transformed_model[i];

		// Uplifting to 4D
		gsl_vector* initial_vec_in_4D = gsl_vector_alloc(4);
		gsl_vector_set(initial_vec_in_4D, 0, gsl_vector_get(initial_vec_in_3D, 0)); // We take X coordinate of 3D vector
		gsl_vector_set(initial_vec_in_4D, 1, gsl_vector_get(initial_vec_in_3D, 1)); // We take Y coordinate of 3D vector
		gsl_vector_set(initial_vec_in_4D, 2, gsl_vector_get(initial_vec_in_3D, 2)); // We take Z coordinate of 3D vector
		gsl_vector_set(initial_vec_in_4D, 3, 1.0); // We use new coordinate for Tau

		// Linear transformation in 4D
		gsl_vector* transformed_vec_in_4D = gsl_vector_alloc(4);
		gsl_blas_dgemv(CblasNoTrans, 1.0, T, initial_vec_in_4D, 0.0, transformed_vec_in_4D);

		// Projection back to 3D
		gsl_vector_set(place_for_transformed_vec_in_3D, 0, gsl_vector_get(transformed_vec_in_4D, 0));
		gsl_vector_set(place_for_transformed_vec_in_3D, 1, gsl_vector_get(transformed_vec_in_4D, 1));
		gsl_vector_set(place_for_transformed_vec_in_3D, 2, gsl_vector_get(transformed_vec_in_4D, 2));

		gsl_vector_free(initial_vec_in_4D);
		gsl_vector_free(transformed_vec_in_4D);
	}
}

void set_model_sequence(Model_3D* model_sequence, gsl_matrix** T_sequence, size_t sequence_length, Model_3D initial_model) {
	for (int j = 0; j < sequence_length; j++) {
		for (int i = 0; i < NUM_OF_VECTORS; i++) {
			model_sequence[j][i] = gsl_vector_alloc(3);
		}
		Model_3D* current_model = &model_sequence[j];
		gsl_matrix* T = T_sequence[j];
		double homotopy_coefficient = (double) j / (double) (sequence_length-1);

		set_transformed_model(*current_model, T, initial_model);
		printf("\nTargeted (x, y, z) coordinates of 3d model points in cm for homotopy coefficient %f:\n", homotopy_coefficient);
		print_model(*current_model);
	}


}
int main(void) {
	const struct NeedleCoordinates state_0 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	const struct NeedleCoordinates state_A = {0.0, 0.0, 0.0, -M_PI/10, M_PI/10, M_PI/10};
	const struct NeedleCoordinates state_B = {0.5, 0.5, 3.0, -M_PI/10, M_PI/10, M_PI/10};

	gsl_matrix* T_0A = gsl_matrix_alloc(4, 4);
	set_transformation_matrix(T_0A, state_0, state_A);

	gsl_matrix* T_AB = gsl_matrix_alloc(4, 4);
	set_transformation_matrix(T_AB, state_A, state_B);

	// Transformation frames are descibed using homotopy Id->T
	// First frame is the the identity transformation (Id).
	// Last frame is the full transformation (T_AB), which brings (state_A) to the (state_B).
	// Other frames are created by mixing T_AB and Id, those frames describe transformation that is "in progress".
	// Example for the case with 3 frames:
	// 	frame 0 is the Id;
	//	frame 1 is 0.5 Id + 0.5 T_AB;
	//	frame 2 is the T_AB;
	const size_t sequence_lenth = 5;
	gsl_matrix* T_Sequence[sequence_lenth];
	set_sequence(T_Sequence, T_AB, sequence_lenth);

	// Array of pointers to vectors that describe 3D model in zero state ({0, 0, 0, 0, 0, 0} state)
	// Model_3D == [p_vec0, p_vec1, p_vec2, ...]
	Model_3D model_in_state_0;
	set_model_in_state_0(model_in_state_0);

	// Array of pointers to vectors that describe 3D model in initial state.
	Model_3D model_in_state_A;
	set_transformed_model(model_in_state_A, T_0A, model_in_state_0);

	// Each 3D model is produced by applying one of transformations "in progress" to the 3D model for initial state.
	// models = [p_model0, p_model1, ...]
	Model_3D model_sequence[sequence_lenth];
	set_model_sequence(model_sequence, T_Sequence, sequence_lenth, model_in_state_A);

 //     gsl_vector* lowest_possible_slider_coordinates[NUM_OF_HEXAGON_VERTICES];
 //     for (int i = 0; i < NUM_OF_HEXAGON_VERTICES; i++) {
 //         lowest_possible_slider_coordinates[i] = gsl_vector_alloc(3);
 //     }
 //     populate_hexagon(lowest_possible_slider_coordinates, BASE_RADIUS, BASE_ANGLE);
 //     printf("\nLowest possible slider coordinates (x, y, z) in cm:\n");
 //     print_lowest_possible_slider_coordinates(lowest_possible_slider_coordinates, NUM_OF_HEXAGON_VERTICES);


/////
//   OptionSliderCoordinates maybe_sliders = find_slider_coordinates(transformed_3dmodel, lowest_possible_slider_coordinates);
//   print_transformed_slider_coordinates(maybe_sliders);

    return 0;
}
