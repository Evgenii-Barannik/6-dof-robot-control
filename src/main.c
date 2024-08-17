#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
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

void print_3dmodel_coordinates(gsl_vector* vectors[NUM_OF_VECTORS]) {
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
    gsl_matrix* T_1to2, // Place where final affine transformation matrix will be written
    const struct NeedleCoordinates state_1,
    const struct NeedleCoordinates state_2
    ) {
	printf("\nInitial coordinates:\n");
	print_needle_coordinates(state_1);
	printf("\nTargeted coordinates:\n");
	print_needle_coordinates(state_2);

	gsl_matrix* T_0to2 = gsl_matrix_alloc(4, 4);
	gsl_matrix_set(T_0to2, 0, 0, cos(state_2.psi) * cos(state_2.theta));
	gsl_matrix_set(T_0to2, 0, 1, cos(state_2.psi) * sin(state_2.theta) * sin(state_2.phi) - sin(state_2.psi) * cos(state_2.phi));
	gsl_matrix_set(T_0to2, 0, 2, cos(state_2.psi) * sin(state_2.theta) * cos(state_2.phi) + sin(state_2.psi) * sin(state_2.phi));
	gsl_matrix_set(T_0to2, 0, 3, state_2.x);

	gsl_matrix_set(T_0to2, 1, 0, sin(state_2.psi) * cos(state_2.theta));
	gsl_matrix_set(T_0to2, 1, 1, sin(state_2.psi) * sin(state_2.theta) * sin(state_2.phi) + cos(state_2.psi) * cos(state_2.phi));
	gsl_matrix_set(T_0to2, 1, 2, sin(state_2.psi) * sin(state_2.theta) * cos(state_2.phi) - cos(state_2.psi) * sin(state_2.phi));
	gsl_matrix_set(T_0to2, 1, 3, state_2.y);

	gsl_matrix_set(T_0to2, 2, 0, -sin(state_2.theta));
	gsl_matrix_set(T_0to2, 2, 1, cos(state_2.theta) * sin(state_2.phi));
	gsl_matrix_set(T_0to2, 2, 2, cos(state_2.theta) * cos(state_2.phi));
	gsl_matrix_set(T_0to2, 2, 3, state_2.z);

	gsl_matrix_set(T_0to2, 3, 0, 0.0);
	gsl_matrix_set(T_0to2, 3, 1, 0.0);
	gsl_matrix_set(T_0to2, 3, 2, 0.0);
	gsl_matrix_set(T_0to2, 3, 3, 1.0);
	printf("\nT_0to2, transformation matrix for (0,0,0,0,0,0) -> (state_2)\n");
	print_matrix_4D(T_0to2);


	gsl_matrix* T_0to1 = gsl_matrix_alloc(4, 4);
	gsl_matrix_set(T_0to1, 0, 0, cos(state_1.psi) * cos(state_1.theta));
	gsl_matrix_set(T_0to1, 0, 1, cos(state_1.psi) * sin(state_1.theta) * sin(state_1.phi) - sin(state_1.psi) * cos(state_1.phi));
	gsl_matrix_set(T_0to1, 0, 2, cos(state_1.psi) * sin(state_1.theta) * cos(state_1.phi) + sin(state_1.psi) * sin(state_1.phi));
	gsl_matrix_set(T_0to1, 0, 3, state_1.x);

	gsl_matrix_set(T_0to1, 1, 0, sin(state_1.psi) * cos(state_1.theta));
	gsl_matrix_set(T_0to1, 1, 1, sin(state_1.psi) * sin(state_1.theta) * sin(state_1.phi) + cos(state_1.psi) * cos(state_1.phi));
	gsl_matrix_set(T_0to1, 1, 2, sin(state_1.psi) * sin(state_1.theta) * cos(state_1.phi) - cos(state_1.psi) * sin(state_1.phi));
	gsl_matrix_set(T_0to1, 1, 3, state_1.y);

	gsl_matrix_set(T_0to1, 2, 0, -sin(state_1.theta));
	gsl_matrix_set(T_0to1, 2, 1, cos(state_1.theta) * sin(state_1.phi));
	gsl_matrix_set(T_0to1, 2, 2, cos(state_1.theta) * cos(state_1.phi));
	gsl_matrix_set(T_0to1, 2, 3, state_1.z);

	gsl_matrix_set(T_0to1, 3, 0, 0.0);
	gsl_matrix_set(T_0to1, 3, 1, 0.0);
	gsl_matrix_set(T_0to1, 3, 2, 0.0);
	gsl_matrix_set(T_0to1, 3, 3, 1.0);
	printf("\nT_0to1, transformation matrix for (0,0,0,0,0,0) -> (state_1)\n");
	print_matrix_4D(T_0to1);


	gsl_matrix* T_0to1_temp = gsl_matrix_alloc(4, 4);
	gsl_matrix_memcpy(T_0to1_temp, T_0to1);

	gsl_permutation * p = gsl_permutation_alloc(4);
	gsl_matrix* T_1to0 = gsl_matrix_alloc(4, 4);
	int signum;
	gsl_linalg_LU_decomp (T_0to1_temp, p, &signum);
	gsl_linalg_LU_invert (T_0to1_temp, p, T_1to0);
	gsl_permutation_free(p);
	printf("\nT_1to0, transformation matrix for (state_1) -> (0,0,0,0,0,0), an inverse of T_0to1\n");
	print_matrix_4D(T_1to0);

	// T_1to2 = T_0to2 * (T_0t1)^(-1)
	// T_1to2 = T_0to2 * T_1to0
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, T_0to2, T_1to0,  0.0, T_1to2);
	printf("\nT_1to2, transformation matrix for (state_1) -> (state_2)\n");
	print_matrix_4D(T_1to2);
}

void set_homotopy_frames(gsl_matrix** homotopy_frames, const gsl_matrix* T, const size_t num_of_homotopy_frames) {
	for (int i = 0; i < num_of_homotopy_frames; i++) {
		double homotopy_coefficient = (double) i / (double) (num_of_homotopy_frames-1);
		homotopy_frames[i] = gsl_matrix_alloc(4, 4);
		gsl_matrix *T_copy = homotopy_frames[i];
		gsl_matrix_memcpy(T_copy, T);

		// Here we mutate affine transformation matrix T into
		// (k)T + (1-k)Id, where k is the homotopy coefficient
		gsl_matrix_scale(T_copy, homotopy_coefficient);
		gsl_matrix* Id = gsl_matrix_alloc(4, 4);
		gsl_matrix_set_identity(Id);
		gsl_matrix_scale(Id, 1.0 - homotopy_coefficient);
		gsl_matrix_add(T_copy, Id);

		printf("\n(%.3f)T_1to2 + (%.3f)Id\n", (homotopy_coefficient), 1-homotopy_coefficient);
		print_matrix_4D(T_copy);
	}
}


void set_initial_3dmodel (gsl_vector** coordinates) {
	for (int i = 0; i < NUM_OF_VECTORS; i++) {
		coordinates[i] = gsl_vector_alloc(3);
	}

	double radius = TABLE_RADIUS;
	double angle_in_radians = TABLE_ANGLE;
	double distance_from_needle_end = NEEDLE_LENGTH;

	gsl_vector_set(coordinates[0], 0, radius * sin(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(coordinates[0], 1, radius * cos(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(coordinates[0], 2, distance_from_needle_end);

	gsl_vector_set(coordinates[1], 0, radius * sin(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(coordinates[1], 1, radius * cos(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(coordinates[1], 2, distance_from_needle_end);

	gsl_vector_set(coordinates[2], 0, radius * sin(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(coordinates[2], 1, radius * cos(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(coordinates[2], 2, distance_from_needle_end);

	gsl_vector_set(coordinates[3], 0, radius * sin(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(coordinates[3], 1, radius * cos(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(coordinates[3], 2, distance_from_needle_end);

	gsl_vector_set(coordinates[4], 0, radius * sin(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(coordinates[4], 1, radius * cos(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(coordinates[4], 2, distance_from_needle_end);

	gsl_vector_set(coordinates[5], 0, radius * sin(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(coordinates[5], 1, radius * cos(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(coordinates[5], 2, distance_from_needle_end);

	// Table centroid
	gsl_vector_set(coordinates[6], 0, 0.0);
	gsl_vector_set(coordinates[6], 1, 0.0);
	gsl_vector_set(coordinates[6], 2, distance_from_needle_end);

	// Needle end
	gsl_vector_set(coordinates[7], 0, 0.0);
	gsl_vector_set(coordinates[7], 1, 0.0);
	gsl_vector_set(coordinates[7], 2, 0.0);

	printf("\nInitial (x, y, z) coordinates of 3d model points in cm:\n");
	print_3dmodel_coordinates(coordinates);
}

// Used to populate slider base coordinates.
// Base coordinates are coordinates of sliders at their lowest possible position at rest.
void set_hexagon(
    gsl_vector** coordinates,
    const double radius,
    const double angle_in_radians) {
    gsl_vector_set(coordinates[0], 0, radius * sin(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[0], 1, radius * cos(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[0], 2, 0.0);

    gsl_vector_set(coordinates[1], 0, radius * sin(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[1], 1, radius * cos(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[1], 2, 0.0);

    gsl_vector_set(coordinates[2], 0, radius * sin(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[2], 1, radius * cos(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[2], 2, 0.0);

    gsl_vector_set(coordinates[3], 0, radius * sin(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[3], 1, radius * cos(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[3], 2, 0.0);

    gsl_vector_set(coordinates[4], 0, radius * sin(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[4], 1, radius * cos(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[4], 2, 0.0);

    gsl_vector_set(coordinates[5], 0, radius * sin(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[5], 1, radius * cos(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[5], 2, 0.0);
}

void set_transformed_coordinates(
    gsl_vector* place_for_transformed_coordinates_in_3D,
    const gsl_matrix* transformation_matrix,
    const gsl_vector* initial_coordinates_in_3D) {

    // Uplifting to 4D
    gsl_vector* initial_coordinates_in_4D = gsl_vector_alloc(4);
    gsl_vector_set(initial_coordinates_in_4D, 0, gsl_vector_get(initial_coordinates_in_3D, 0)); // We take X coordinate of 3D vector
    gsl_vector_set(initial_coordinates_in_4D, 1, gsl_vector_get(initial_coordinates_in_3D, 1)); // We take Y coordinate of 3D vector
    gsl_vector_set(initial_coordinates_in_4D, 2, gsl_vector_get(initial_coordinates_in_3D, 2)); // We take Z coordinate of 3D vector
    gsl_vector_set(initial_coordinates_in_4D, 3, 1.0); // We use new coordinate for Tau

    // Linear transformation in 4D
    gsl_vector* transformed_coordinates_in_4D = gsl_vector_alloc(4);
    gsl_blas_dgemv(CblasNoTrans, 1.0, transformation_matrix, initial_coordinates_in_4D, 0.0, transformed_coordinates_in_4D);

    // Projection back to 3D
    gsl_vector_set(place_for_transformed_coordinates_in_3D, 0, gsl_vector_get(transformed_coordinates_in_4D, 0));
    gsl_vector_set(place_for_transformed_coordinates_in_3D, 1, gsl_vector_get(transformed_coordinates_in_4D, 1));
    gsl_vector_set(place_for_transformed_coordinates_in_3D, 2, gsl_vector_get(transformed_coordinates_in_4D, 2));

    gsl_vector_free(initial_coordinates_in_4D);
    gsl_vector_free(transformed_coordinates_in_4D);
}

OptionSliderCoordinates find_slider_coordinates(
    gsl_vector** table_transformed_coordinates,
    gsl_vector** base_coordinates) {

    OptionSliderCoordinates maybe_slider_coordinates;

    for (int i = 0; i < NUM_OF_SLIDERS; i++) {
        maybe_slider_coordinates.coordinates[i] = gsl_vector_alloc(3);
        gsl_vector_memcpy(maybe_slider_coordinates.coordinates[i], base_coordinates[i]);

        double table_transformed_x = gsl_vector_get(table_transformed_coordinates[i], 0);
        double table_transformed_y = gsl_vector_get(table_transformed_coordinates[i], 1);
        double table_transformed_z = gsl_vector_get(table_transformed_coordinates[i], 2);
        double base_x = gsl_vector_get(base_coordinates[i], 0);
        double base_y = gsl_vector_get(base_coordinates[i], 1);

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
                gsl_vector_set(maybe_slider_coordinates.coordinates[i], 2, slider_z);
            } else {
                gsl_vector_free(maybe_slider_coordinates.coordinates[i]);
                maybe_slider_coordinates.coordinates[i] = NULL;
            };
        } else {
            gsl_vector_free(maybe_slider_coordinates.coordinates[i]);
            maybe_slider_coordinates.coordinates[i] = NULL;
        }
    }
    return maybe_slider_coordinates;
}

typedef gsl_matrix HomotopyFrame;
typedef gsl_vector* Model[NUM_OF_VECTORS];

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
	assert(velocity_max >= 0.0);
	assert(progress >= 0.0 && progress <= 1.0);
	assert(ACCURACY >= 0.0 && ACCURACY <= 1.0);

	double k = log((1.0-ACCURACY)/(ACCURACY));
	double velocity = velocity_max/(1.0 + exp(k*(2.0*progress - 1.0)));
	return velocity;
}

double targeted_acceleration(const double progress, const double velocity_max) {
	// https://www.desmos.com/calculator/j0sgzmoqbj
	assert(velocity_max >= 0.0);
	assert(progress >= 0.0 && progress <= 1.0);
	assert(ACCURACY >= 0.0 && ACCURACY <= 1.0);

	double k = log((1.0-ACCURACY)/(ACCURACY));
	double acceleration = -2*k*velocity_max * exp(k*(2.0*progress - 1.0)) / pow(exp(k*(2.0*progress - 1.0))+1, 2);
	return acceleration;
}

int main(void) {
	printf("%f", targeted_acceleration(0.5, 1.0));

	const struct NeedleCoordinates initial_state = {0.0, 0.0, 0.0, -M_PI/10, M_PI/10, M_PI/10};
	const struct NeedleCoordinates targeted_state = {0.5, 0.5, 3.0, -M_PI/10, M_PI/10, M_PI/10};
	const size_t num_of_homotopy_frames = 5;
	gsl_matrix* T = gsl_matrix_alloc(4, 4);
	set_transformation_matrix(T, initial_state, targeted_state);

	// Frames for homotopy Id->T
	// First frame is the the identity transformation (Id).
	// Last frame is the full transformation (T), it brings (initial_state) to the (targeted_state).
	// Other frames are created by mixing T and Id, they describe affine transformation that is "in progress".
	// Example for case with 3 frames in total:
	// 	frame 0 is the Id;
	//	frame 1 is 0.5 Id + 0.5 T;
	//	frame 2 is the T;
	HomotopyFrame* homotopy_frames[num_of_homotopy_frames];
	set_homotopy_frames(homotopy_frames, T, num_of_homotopy_frames);

	// Array of pointers to vectors that describe initial 3D model.
	// Model == [p_vec0, p_vec1, p_vec2, ...]
	Model initial_3dmodel;
	set_initial_3dmodel(initial_3dmodel);

	// Each 3D modes is produced by applying one of transformations "in progress".
	// models = [p_model0, p_model1, ...]
	Model* models[num_of_homotopy_frames];
	for (int j = 0; j < num_of_homotopy_frames; j++) {
		models[j] = malloc(sizeof(Model));
		Model* model = models[j];
		HomotopyFrame* frame = homotopy_frames[j];
		double homotopy_coefficient = (double) j / (double) (num_of_homotopy_frames-1);
		printf("\nTargeted (x, y, z) coordinates of 3d model points in cm for homotopy coefficient %f:\n", homotopy_coefficient);
			
		for (int i = 0; i < NUM_OF_VECTORS; i++) {
			(*model)[i] = gsl_vector_alloc(3);
			const gsl_vector* initial_vec= initial_3dmodel[i];
			gsl_vector* vec = (*model)[i];
			set_transformed_coordinates(vec, frame, initial_vec);
		}
		print_3dmodel_coordinates(*model);
	}

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
