#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_blas.h>

const double DEGREES_TO_RADIANS = (2 * M_PI / 360);
const int NUM_OF_HEXAGON_VERTICES = 6;

const double TABLE_RADIUS = 3.0; // Radius of the hexagonal table in cm.
const double TABLE_ANGLE = DEGREES_TO_RADIANS * 80; // Angle between vertex pairs of the hexagonal table in radians. A regular hexagon will have 60 deg.
const double BASE_RADIUS = 7.0; // Radius of the hexagonal base in cm.
const double BASE_ANGLE = DEGREES_TO_RADIANS * 15; // Angle between vertex pairs of the hexagonal base in radians. A regular hexagon will have 60 deg.
const double MAX_POSSIBLE_HEIGHT = 15.0; // Highest possible position for sliders in cm.
const double MIN_POSSIBLE_HEIGHT = 0; // Lowest possible position for sliders in cm.
const double ARM_LENGTH = 8.0; // Length of all robot arms in cm.

struct TableCoordinates { // Full description of the table position.
    double x, y, z; // Coordinates of the table center.
    double phi, theta, psi; // Euler angles for the table.
};

// Array of pointers to NULL (if position for any slider is not possible)
// or pointers to slider positions (if all slider positions are possible).
typedef struct {
    gsl_vector* positions[NUM_OF_HEXAGON_VERTICES];
} OptionSliderPositions;

bool all_not_NULL(const OptionSliderPositions* options) {
    for (int i = 0; i < NUM_OF_HEXAGON_VERTICES; ++i) {
        gsl_vector* vec_p = options->positions[i];
        if (vec_p == NULL) {
            return false;
        }
    }
    return true;
}

void print_TableCoordinates(struct TableCoordinates table) {
    printf("\n");
    printf("Targeted table coordinates: \n");
    printf("Centroid (x, y, z) coordinates in cm:      (%f, %f, %f)\n", table.x, table.y, table.z);
    printf("Euler angles (phi, theta, psi) in radians: (%f, %f, %f)\n", table.phi, table.theta, table.psi);
}

void print_points(gsl_vector** points, int num_points) {
    if (num_points == NUM_OF_HEXAGON_VERTICES || num_points == NUM_OF_HEXAGON_VERTICES + 1) {
            for(int i = 0; i < num_points; i++) {
                printf("%d: (%f, %f, %f)", i,
                    gsl_vector_get(points[i], 0),
                    gsl_vector_get(points[i], 1),
                    gsl_vector_get(points[i], 2));

                if (i==NUM_OF_HEXAGON_VERTICES) {
                    printf(" <-- centroid");
                }
                printf("\n");
            }
    } else {
        printf("Error, number of points is incorrect");
    }
}


void print_slider_positions(const OptionSliderPositions option) {
    bool all_positions_found = all_not_NULL(&option);

    if (all_positions_found) {
        printf("\nTargeted table coordinates can be realized using these slider positions:\n");
        for (int i = 0; i < NUM_OF_HEXAGON_VERTICES; i++) {
            printf("%d: (%f, %f, %f)\n", i,
                gsl_vector_get(option.positions[i], 0),
                gsl_vector_get(option.positions[i], 1),
                gsl_vector_get(option.positions[i], 2));
        }
    } else {
        printf("\nFailed to find slider positions that can realize targeted table coordinates.\n");
        printf("Some SliderPositions are pointing to NULL:\n");
        for (int i = 0; i < NUM_OF_HEXAGON_VERTICES; i++) {
            bool if_current_position_found = (option.positions[i] != NULL);
            if (if_current_position_found) {
                printf("%d: ok \n", i);
            } else {
                printf("%d: NULL\n", i);
            }
        }
    }
    printf("\n");
}

void populate_transformation_matrix(
    gsl_matrix * transformation_matrix,
    const struct TableCoordinates c) {

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

void populate_hexagon(
    gsl_vector **points,
    const double radius,
    const double angle_in_radians) {

    gsl_vector_set(points[0], 0, radius * sin(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
    gsl_vector_set(points[0], 1, radius * cos(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
    gsl_vector_set(points[0], 2, 0.0);

    gsl_vector_set(points[1], 0, radius * sin(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
    gsl_vector_set(points[1], 1, radius * cos(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
    gsl_vector_set(points[1], 2, 0.0);

    gsl_vector_set(points[2], 0, radius * sin(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
    gsl_vector_set(points[2], 1, radius * cos(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
    gsl_vector_set(points[2], 2, 0.0);

    gsl_vector_set(points[3], 0, radius * sin(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
    gsl_vector_set(points[3], 1, radius * cos(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
    gsl_vector_set(points[3], 2, 0.0);

    gsl_vector_set(points[4], 0, radius * sin(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
    gsl_vector_set(points[4], 1, radius * cos(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
    gsl_vector_set(points[4], 2, 0.0);

    gsl_vector_set(points[5], 0, radius * sin(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
    gsl_vector_set(points[5], 1, radius * cos(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
    gsl_vector_set(points[5], 2, 0.0);
}

void transform_point(
    gsl_vector *for_transformed_point_in_3D,
    const gsl_matrix *transformation_matrix,
    const gsl_vector *initial_point_in_3D) {

    // Uplifting to 4D
    gsl_vector *initial_point_in_4D = gsl_vector_alloc(4);
    gsl_vector_set(initial_point_in_4D, 0, gsl_vector_get(initial_point_in_3D, 0));
    gsl_vector_set(initial_point_in_4D, 1, gsl_vector_get(initial_point_in_3D, 1));
    gsl_vector_set(initial_point_in_4D, 2, gsl_vector_get(initial_point_in_3D, 2));
    gsl_vector_set(initial_point_in_4D, 3, 1.0);

    // Linear transformation in 4D
    gsl_vector *transformed_point_in_4D = gsl_vector_alloc(4);
    gsl_blas_dgemv(CblasNoTrans, 1.0, transformation_matrix, initial_point_in_4D, 0.0, transformed_point_in_4D);

    // Projection back to 3D
    gsl_vector_set(for_transformed_point_in_3D, 0, gsl_vector_get(transformed_point_in_4D, 0));
    gsl_vector_set(for_transformed_point_in_3D, 1, gsl_vector_get(transformed_point_in_4D, 1));
    gsl_vector_set(for_transformed_point_in_3D, 2, gsl_vector_get(transformed_point_in_4D, 2));

    gsl_vector_free(initial_point_in_4D);
    gsl_vector_free(transformed_point_in_4D);
}

void populate_centroid(gsl_vector **points) {
    double x_sum = 0.0, y_sum = 0.0, z_sum = 0.0;
    for (int j = 0; j < NUM_OF_HEXAGON_VERTICES; j++) {
        x_sum += gsl_vector_get(points[j], 0);
        y_sum += gsl_vector_get(points[j], 1);
        z_sum += gsl_vector_get(points[j], 2);

    }
    gsl_vector_set(points[NUM_OF_HEXAGON_VERTICES], 0, x_sum/NUM_OF_HEXAGON_VERTICES);
    gsl_vector_set(points[NUM_OF_HEXAGON_VERTICES], 1, y_sum/NUM_OF_HEXAGON_VERTICES);
    gsl_vector_set(points[NUM_OF_HEXAGON_VERTICES], 2, z_sum/NUM_OF_HEXAGON_VERTICES);
}

OptionSliderPositions find_slider_positions(
    gsl_vector **table_transformed_points,
    gsl_vector **base_points) {

    OptionSliderPositions maybe_slider_positions;

    for (int i = 0; i < NUM_OF_HEXAGON_VERTICES; i++) {
        maybe_slider_positions.positions[i] = gsl_vector_alloc(3);
        gsl_vector_memcpy(maybe_slider_positions.positions[i], base_points[i]);

        double table_transformed_point_x = gsl_vector_get(table_transformed_points[i], 0);
        double table_transformed_point_y = gsl_vector_get(table_transformed_points[i], 1);
        double table_transformed_point_z = gsl_vector_get(table_transformed_points[i], 2);
        double base_point_x = gsl_vector_get(base_points[i], 0);
        double base_point_y = gsl_vector_get(base_points[i], 1);

        bool slider_inside_range = false;
        double z_part = pow(ARM_LENGTH, 2) - pow((base_point_x - table_transformed_point_x), 2) - pow((base_point_y - table_transformed_point_y), 2);

        if (z_part > 0) {
            // This is one of two possible solutions for z coordinate of slider (higher solution of two).
            // Lower solution will be table_transformed_point_z - sqrt(z_part).
            // We don't use second solution because we fear damage during solution switching.
            double z_coordinate_of_slider = table_transformed_point_z + sqrt(z_part);
            slider_inside_range = (z_coordinate_of_slider < MAX_POSSIBLE_HEIGHT) && (z_coordinate_of_slider > MIN_POSSIBLE_HEIGHT);
            if (slider_inside_range) {
                gsl_vector_set(maybe_slider_positions.positions[i], 2, z_coordinate_of_slider);
            } else {
                gsl_vector_free(maybe_slider_positions.positions[i]);
                maybe_slider_positions.positions[i] = NULL;
            };
        } else {
            gsl_vector_free(maybe_slider_positions.positions[i]);
            maybe_slider_positions.positions[i] = NULL;
        }
    }
    return maybe_slider_positions;
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
    struct TableCoordinates targeted_table_coords = {0.5, -0.5, 5, -M_PI/10, M_PI/3, 0};

    gsl_matrix *transformation_matrix = gsl_matrix_alloc(4, 4);
    populate_transformation_matrix(transformation_matrix, targeted_table_coords);

    // Array with vertices plus an empty place for centroid:
    gsl_vector* initial_table_points[NUM_OF_HEXAGON_VERTICES+1];
    for (int i = 0; i < NUM_OF_HEXAGON_VERTICES+1; i++) {
        initial_table_points[i] = gsl_vector_alloc(3);
    }
    populate_hexagon(initial_table_points, TABLE_RADIUS, TABLE_ANGLE); // Populates all vectors except last
    populate_centroid(initial_table_points); // Populates last vector

    // Transformation of all points (including centroid):
    gsl_vector *transformed_table_points[NUM_OF_HEXAGON_VERTICES+1];
    for (int i = 0; i < NUM_OF_HEXAGON_VERTICES+1; i++) {
        transformed_table_points[i] = gsl_vector_alloc(3);
        transform_point(transformed_table_points[i], transformation_matrix, initial_table_points[i]);
    }

    printf("\nInitial (x, y, z) coordinates of table points in cm, with vertices and centroid:\n");
    print_points(initial_table_points, NUM_OF_HEXAGON_VERTICES+1);

    printf("\nTargeted (x, y, z) coordinates of table points in cm, with vertices and centroid:\n");
    print_points(transformed_table_points, NUM_OF_HEXAGON_VERTICES+1);

    gsl_vector *base_points[NUM_OF_HEXAGON_VERTICES];
    for (int i = 0; i < NUM_OF_HEXAGON_VERTICES; i++) {
        base_points[i] = gsl_vector_alloc(3);
    }
    populate_hexagon(base_points, BASE_RADIUS, BASE_ANGLE); // Each base point lay exactly under corresponing slider.

    printf("\nImmutable (x, y, z) coordinates of base points in cm, with vertices only:\n");
    print_points(base_points, NUM_OF_HEXAGON_VERTICES);

    print_TableCoordinates(targeted_table_coords);

    OptionSliderPositions maybe_sliders = find_slider_positions(transformed_table_points, base_points);
    print_slider_positions(maybe_sliders);

    return 0;
}
