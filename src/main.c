#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_blas.h>

// Ensures that assertions are not disabled
#undef NDEBUG
#include <assert.h>


/// CONSTANTS THAT SHOULD NOT BE CHANGED
const double DEGREES_TO_RADIANS = (2 * M_PI / 360);
const int NUM_OF_HEXAGON_VERTICES = 6;
const int NUM_OF_TRANSFORMED_VECTORS = NUM_OF_HEXAGON_VERTICES + 2;
const int NUM_OF_SLIDERS = 6;
///


/// CONSTANTS THAT CAN BE CHANGED
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
    gsl_vector* coordinates[NUM_OF_HEXAGON_VERTICES];
} OptionSliderCoordinates;

bool all_not_NULL(const OptionSliderCoordinates* options) {
    for (int i = 0; i < NUM_OF_HEXAGON_VERTICES; ++i) {
        gsl_vector* vec_p = options->coordinates[i];
        if (vec_p == NULL) {
            return false;
        }
    }
    return true;
}

void print_targeted_needle_coordinates(struct NeedleCoordinates table) {
    printf("\n");
    printf("Targeted needle coordinates: \n");
    printf("Linear coordinates (x, y, z) in cm:        (%f, %f, %f)\n", table.x, table.y, table.z);
    printf("Euler angles (phi, theta, psi) in radians: (%f, %f, %f)\n", table.phi, table.theta, table.psi);
}

void print_3dmodel_coordinates(gsl_vector** vectors, int num_of_vectors_to_print) {
    assert(num_of_vectors_to_print == NUM_OF_HEXAGON_VERTICES + 2);

    for(int i = 0; i < num_of_vectors_to_print; i++) {
        printf("%d: (%f, %f, %f)", i,
            gsl_vector_get(vectors[i], 0),
            gsl_vector_get(vectors[i], 1),
            gsl_vector_get(vectors[i], 2));

        if (i == NUM_OF_HEXAGON_VERTICES) {
            printf(" <-- centroid");
        }
        if (i == NUM_OF_HEXAGON_VERTICES + 1){
            printf(" <-- needle");
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
        for (int i = 0; i < NUM_OF_HEXAGON_VERTICES; i++) {
            printf("Slider %d: (%f, %f, %f)\n", i,
                gsl_vector_get(option.coordinates[i], 0),
                gsl_vector_get(option.coordinates[i], 1),
                gsl_vector_get(option.coordinates[i], 2));
        }
    } else {
        printf("\n=== Solution NOT found! ===\n");
        printf("Failed to find slider coordinates that can realize targeted table coordinates.\n");
        printf("Targeted table coordinates deemed infeasible. Slider diagnostics:\n");
        for (int i = 0; i < NUM_OF_HEXAGON_VERTICES; i++) {
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

void populate_transformation_matrix(
    gsl_matrix* transformation_matrix,
    const struct NeedleCoordinates c) {
    gsl_matrix_set(transformation_matrix, 0, 0, cos(c.psi) * cos(c.theta));
    gsl_matrix_set(transformation_matrix, 0, 1, cos(c.psi) * sin(c.theta) * sin(c.phi) - sin(c.psi) * cos(c.phi));
    gsl_matrix_set(transformation_matrix, 0, 2, cos(c.psi) * sin(c.theta) * cos(c.phi) + sin(c.psi) * sin(c.phi));
    gsl_matrix_set(transformation_matrix, 0, 3, c.x);

    gsl_matrix_set(transformation_matrix, 1, 0, sin(c.psi) * cos(c.theta));
    gsl_matrix_set(transformation_matrix, 1, 1, sin(c.psi) * sin(c.theta) * sin(c.phi) + cos(c.psi) * cos(c.phi));
    gsl_matrix_set(transformation_matrix, 1, 2, sin(c.psi) * sin(c.theta) * cos(c.phi) - cos(c.psi) * sin(c.phi));
    gsl_matrix_set(transformation_matrix, 1, 3, c.y);

    gsl_matrix_set(transformation_matrix, 2, 0, -sin(c.theta));
    gsl_matrix_set(transformation_matrix, 2, 1, cos(c.theta) * sin(c.phi));
    gsl_matrix_set(transformation_matrix, 2, 2, cos(c.theta) * cos(c.phi));
    gsl_matrix_set(transformation_matrix, 2, 3, c.z);

    gsl_matrix_set(transformation_matrix, 3, 0, 0.0);
    gsl_matrix_set(transformation_matrix, 3, 1, 0.0);
    gsl_matrix_set(transformation_matrix, 3, 2, 0.0);
    gsl_matrix_set(transformation_matrix, 3, 3, 1.0);
}

void populate_table_3dmodel (gsl_vector** coordinates) {
    double radius = TABLE_RADIUS;
    double angle_in_radians = TABLE_ANGLE;
    double separation_from_needle_end = NEEDLE_LENGTH;

    gsl_vector_set(coordinates[0], 0, radius * sin(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[0], 1, radius * cos(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[0], 2, separation_from_needle_end);

    gsl_vector_set(coordinates[1], 0, radius * sin(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[1], 1, radius * cos(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[1], 2, separation_from_needle_end);

    gsl_vector_set(coordinates[2], 0, radius * sin(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[2], 1, radius * cos(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[2], 2, separation_from_needle_end);

    gsl_vector_set(coordinates[3], 0, radius * sin(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[3], 1, radius * cos(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[3], 2, separation_from_needle_end);

    gsl_vector_set(coordinates[4], 0, radius * sin(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[4], 1, radius * cos(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[4], 2, separation_from_needle_end);

    gsl_vector_set(coordinates[5], 0, radius * sin(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[5], 1, radius * cos(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
    gsl_vector_set(coordinates[5], 2, separation_from_needle_end);

    // Table centroid
    gsl_vector_set(coordinates[6], 0, 0.0);
    gsl_vector_set(coordinates[6], 1, 0.0);
    gsl_vector_set(coordinates[6], 2, separation_from_needle_end);

    // Needle end
    gsl_vector_set(coordinates[7], 0, 0.0);
    gsl_vector_set(coordinates[7], 1, 0.0);
    gsl_vector_set(coordinates[7], 2, 0.0);
}

// Used to populate both table and base coordinates.
// Base coordinates are coordinates of sliders at their lowest possible position at rest.
void populate_hexagon(
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

void populate_transformed_coordinates(
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

    for (int i = 0; i < NUM_OF_HEXAGON_VERTICES; i++) {
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
// Below is an example of how damage during solution switching may occur.
// Damage during solution switching appears unlikely when considering only one solution.
// However, it becomes possible when both solutions are calculated simultaneously.
// Lets imagine that for given targeted table coordinates there are
// two possible positions (lower and higher ones) for each of our 6 sliders.
// Lets imagine that each slider is present in its higher possible position.
// Lets now imagine that for new targeted table coordinates our slider #3
// has only one possible solution left (lower one).
// So the higher solution is now impossible, but the lower one is still possible.
// If we are not forcing each slider step to be infinitely small,
// our system may decide that slider #3 should now switch to the lower position,
// thus passing through set of incorrect configurations.
// This movement can damage our mechanical system.

int main(void) {
    struct NeedleCoordinates targeted_needle_coordinates = {0.5, 0.5, 3.0, -M_PI/10, M_PI/10, M_PI/10};

    gsl_matrix* transformation_matrix = gsl_matrix_alloc(4, 4);
    // Transformation that brings point (0, 0, 0, 0, 0 ,0) NeedleCoordinates to the targeted_needle_coordinates.
    populate_transformation_matrix(transformation_matrix, targeted_needle_coordinates);

    // Array of all vectors describing our 3D model.
    gsl_vector* initial_3dmodel[NUM_OF_TRANSFORMED_VECTORS];
    for (int i = 0; i < NUM_OF_TRANSFORMED_VECTORS; i++) {
        initial_3dmodel[i] = gsl_vector_alloc(3);
    }
    populate_table_3dmodel(initial_3dmodel);

    // Transformation of all vectors:
    gsl_vector* transformed_3dmodel[NUM_OF_TRANSFORMED_VECTORS];
    for (int i = 0; i < NUM_OF_TRANSFORMED_VECTORS; i++) {
        transformed_3dmodel[i] = gsl_vector_alloc(3);
        populate_transformed_coordinates(transformed_3dmodel[i], transformation_matrix, initial_3dmodel[i]);
    }

    print_targeted_needle_coordinates(targeted_needle_coordinates);

    printf("\nInitial (x, y, z) coordinates of the table 3d model points in cm:\n");
    print_3dmodel_coordinates(initial_3dmodel, NUM_OF_TRANSFORMED_VECTORS);

    printf("\nTargeted (x, y, z) coordinates of the table 3d model points in cm:\n");
    print_3dmodel_coordinates(transformed_3dmodel, NUM_OF_TRANSFORMED_VECTORS);

    gsl_vector* lowest_possible_slider_coordinates[NUM_OF_HEXAGON_VERTICES];
    for (int i = 0; i < NUM_OF_HEXAGON_VERTICES; i++) {
        lowest_possible_slider_coordinates[i] = gsl_vector_alloc(3);
    }
    populate_hexagon(lowest_possible_slider_coordinates, BASE_RADIUS, BASE_ANGLE);

    printf("\nLowest possible slider coordinates (x, y, z) in cm:\n");
    print_lowest_possible_slider_coordinates(lowest_possible_slider_coordinates, NUM_OF_HEXAGON_VERTICES);

    OptionSliderCoordinates maybe_sliders = find_slider_coordinates(transformed_3dmodel, lowest_possible_slider_coordinates);
    print_transformed_slider_coordinates(maybe_sliders);

    return 0;
}
