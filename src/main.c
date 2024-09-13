#include <stddef.h>
#include <string.h> // for memcpy()
#include <stdbool.h> // for bool
#include <time.h>   // for clock_gettime()
#include <unistd.h> // for usleep()
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#undef NDEBUG // Ensures that assertions are not disabled
#include <assert.h> // for assert()


/// CONSTANTS THAT SHOULD NOT BE CHANGED
const double DEGREES_TO_RADIANS = (2 * M_PI / 360);
const size_t NUM_OF_HEXAGON_VERTICES = 6;
const size_t NUM_OF_VECTORS = NUM_OF_HEXAGON_VERTICES + 2;
const size_t NUM_OF_SLIDERS = 6;
///


/// CONSTANTS THAT CAN BE CHANGED
const bool PRINT_DEBUG_INFO_GLOBAL = true;
const size_t NUM_OF_FRAMES_TO_PRINT_ON_CREATION = 2;  // Must be >= 1
const size_t NUM_OF_FRAMES_TO_PRINT_ON_FOLLOWING = 2; // Must be >= 1
const size_t NUM_OF_FRAMES = 10 * 1000; // Must be >= 2
const size_t MOVEMENT_DURATION_MKS = 5 * 1000 * 1000;
const double MIN_SLIDER_STEP_CM = 0.001;
const double NEEDLE_LENGTH_CM = 3.0;
const double TABLE_RADIUS_CM = 2.5; // Radius of the hexagonal table in cm.
const double SLIDERS_RADIUS_CM = 7.0; // Radius of the slider hexagonal in cm.
// Angle between vertices (0, 1) and (2, 3) and (4, 5) of the hexagonal table in radians. A regular hexagon will have 60 deg:
const double TABLE_ANGLE_RAD = DEGREES_TO_RADIANS * 80; 
// Angle between sliders (0, 1) and (2, 3) and (4, 5) in radians. Look at any horizontal plane to see this angle. A regular hexagon will have 60 deg:
const double SLIDERS_ANGLE_RAD = DEGREES_TO_RADIANS * 15; 
const double MAX_SLIDER_HEIGHT_CM = 20.0; // Highest possible position for sliders in cm.
const double MIN_SLIDER_HEIGHT_CM = 0; // Lowest possible position for sliders in cm.
const double SLIDER_ARM_LENGTH_CM = 8.0; // Length of all robot arms in cm.
///


// Full description of the needle position, and thus, whole system ideally.
typedef struct { 
	double x, y, z; // Coordinates of the needle end in cm.
	double roll, pitch, yaw; // Euler angles for the needle in radians.
} NeedleCoordinates;

typedef struct {
    gsl_vector* vecs [NUM_OF_VECTORS];
} Table_Model;

// Array of pointers to NULL (if sliders were not found)
// or pointers to slider vectors (if sliders were found).
typedef struct {
    gsl_vector* vecs [NUM_OF_SLIDERS];
} Sliders_Model_Nullable;

typedef struct {
	bool slider_movement_needed [NUM_OF_SLIDERS];
	bool slider_movement_direction_is_up [NUM_OF_SLIDERS];
} Sliders_Movement_Recipe;

typedef struct { 
	Table_Model table_frames [NUM_OF_FRAMES]; 
	Sliders_Model_Nullable sliders_frames [NUM_OF_FRAMES];
	uint64_t timestamps [NUM_OF_FRAMES];
	size_t num_of_frames;

	// Pointers to single structs
	NeedleCoordinates* initial_state;
	NeedleCoordinates* final_state;
} Trajectory;

void print_needle_coordinates(const NeedleCoordinates* state) {
	printf("Linear coordinates (x, y, z) in cm:      (%10.6f, %10.6f, %10.6f)\n", state->x, state->y, state->z);
	printf("Euler angles (roll, pitch, yaw) in rads: (%10.6f, %10.6f, %10.6f)\n", state->roll, state->pitch, state->yaw);
}

void print_table(const Table_Model* table) {
	for(size_t i = 0; i < NUM_OF_VECTORS; i++) {

		printf("%lu: (%10.6f, %10.6f, %10.6f)", i,
				gsl_vector_get(table->vecs[i], 0),
				gsl_vector_get(table->vecs[i], 1),
				gsl_vector_get(table->vecs[i], 2));

		if (i == NUM_OF_VECTORS - 2) {
			printf(" <-- table centroid");
		}
		if (i == NUM_OF_VECTORS - 1){
			printf(" <-- needle end");
		}
		printf("\n");
	}
}

void print_sliders(const Sliders_Model_Nullable* sliders) {
	for(size_t i = 0; i < NUM_OF_SLIDERS; i++) {
		if (sliders->vecs[i] != NULL) {
			printf("Slider %lu: (%10.6f, %10.6f, %10.6f)\n", i,
					gsl_vector_get(sliders->vecs[i], 0),
					gsl_vector_get(sliders->vecs[i], 1),
					gsl_vector_get(sliders->vecs[i], 2));
		} else {
			printf("Slider %lu: ERROR\n", i);
		}
	}
}

Sliders_Model_Nullable get_slider_difference(const Sliders_Model_Nullable sliders_2, const Sliders_Model_Nullable sliders_1) {
    Sliders_Model_Nullable difference;

    for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
        if (sliders_1.vecs[i] != NULL && sliders_2.vecs[i] != NULL) {
            difference.vecs[i] = gsl_vector_alloc(3);
            
            for (size_t j = 0; j < 3; j++) {
                double diff = gsl_vector_get(sliders_2.vecs[i], j) - gsl_vector_get(sliders_1.vecs[i], j);
                gsl_vector_set(difference.vecs[i], j, diff);
            }
        } else {
            difference.vecs[i] = NULL;
        }
    }

    return difference;
}

void print_trajectory(const Trajectory* trajectory, const bool print_frames) {
	const size_t num_of_frames = trajectory->num_of_frames;
	printf("\n==== Trajectory information:");
	printf("\nNumber of frames: %lu", num_of_frames);
	printf("\nTrajectory start timestamp: %llu microseconds\n", trajectory->timestamps[0]);
	printf("Trajectory end timestamp:   %llu microseconds\n", trajectory->timestamps[num_of_frames - 1]);


	printf("\n==== Needle coordinates for the initial_state:\n");
	print_needle_coordinates(trajectory->initial_state);
	printf("\n==== Needle coordinates for the final_state:\n");
	print_needle_coordinates(trajectory->final_state);

	if (print_frames) {
		size_t print_interval = (num_of_frames + NUM_OF_FRAMES_TO_PRINT_ON_CREATION - 1) / NUM_OF_FRAMES_TO_PRINT_ON_CREATION;
		for (size_t i = 0; i < num_of_frames; i+= print_interval) {
			double progress = (double) i / (double) (num_of_frames - 1);
			printf("\n==== Frame %lu / %lu\n", i, num_of_frames - 1);
			printf("Trajectory progress: %.3f%%\n", progress * 100.0);
			printf("Timestamp: %llu microseconds\n", trajectory->timestamps[i]);
			printf("\nTable coordinates (x, y, z) in cm:\n");
			print_table(&trajectory->table_frames[i]);
			printf("\nSlider coordinates (x, y, z) in cm:\n");
			print_sliders(&trajectory->sliders_frames[i]);
		}
	}
}


void print_sliders_movement_recipe(const Sliders_Movement_Recipe* recipe) {
    printf("\nSliders Movement Recipe:\n");
	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
		printf("Slider %lu movement needed: %s\n", i, recipe->slider_movement_needed[i] ? "True" : "False");

	}
	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
		printf("Slider %lu movement direction is upwards: %s\n", i, recipe->slider_movement_direction_is_up[i] ? "True" : "False");
	}
}

gsl_matrix* get_transformation_matrix(const NeedleCoordinates state_A, const NeedleCoordinates state_B) {
	// M_TB, transformation matrix for (state_top) -> (state_B)
	gsl_matrix* M_TB = gsl_matrix_alloc(4, 4);
	gsl_matrix_set(M_TB, 0, 0, cos(state_B.yaw) * cos(state_B.pitch));
	gsl_matrix_set(M_TB, 0, 1, cos(state_B.yaw) * sin(state_B.pitch) * sin(state_B.roll) - sin(state_B.yaw) * cos(state_B.roll));
	gsl_matrix_set(M_TB, 0, 2, cos(state_B.yaw) * sin(state_B.pitch) * cos(state_B.roll) + sin(state_B.yaw) * sin(state_B.roll));
	gsl_matrix_set(M_TB, 0, 3, state_B.x);
	gsl_matrix_set(M_TB, 1, 0, sin(state_B.yaw) * cos(state_B.pitch));
	gsl_matrix_set(M_TB, 1, 1, sin(state_B.yaw) * sin(state_B.pitch) * sin(state_B.roll) + cos(state_B.yaw) * cos(state_B.roll));
	gsl_matrix_set(M_TB, 1, 2, sin(state_B.yaw) * sin(state_B.pitch) * cos(state_B.roll) - cos(state_B.yaw) * sin(state_B.roll));
	gsl_matrix_set(M_TB, 1, 3, state_B.y);
	gsl_matrix_set(M_TB, 2, 0, -sin(state_B.pitch));
	gsl_matrix_set(M_TB, 2, 1, cos(state_B.pitch) * sin(state_B.roll));
	gsl_matrix_set(M_TB, 2, 2, cos(state_B.pitch) * cos(state_B.roll));
	gsl_matrix_set(M_TB, 2, 3, state_B.z);
	gsl_matrix_set(M_TB, 3, 0, 0.0);
	gsl_matrix_set(M_TB, 3, 1, 0.0);
	gsl_matrix_set(M_TB, 3, 2, 0.0);
	gsl_matrix_set(M_TB, 3, 3, 1.0);

	// M_TA, transformation matrix for (state_top) -> (state_A)
	gsl_matrix* M_TA = gsl_matrix_alloc(4, 4);
	gsl_matrix_set(M_TA, 0, 0, cos(state_A.yaw) * cos(state_A.pitch));
	gsl_matrix_set(M_TA, 0, 1, cos(state_A.yaw) * sin(state_A.pitch) * sin(state_A.roll) - sin(state_A.yaw) * cos(state_A.roll));
	gsl_matrix_set(M_TA, 0, 2, cos(state_A.yaw) * sin(state_A.pitch) * cos(state_A.roll) + sin(state_A.yaw) * sin(state_A.roll));
	gsl_matrix_set(M_TA, 0, 3, state_A.x);
	gsl_matrix_set(M_TA, 0, 3, state_A.x);
	gsl_matrix_set(M_TA, 1, 0, sin(state_A.yaw) * cos(state_A.pitch));
	gsl_matrix_set(M_TA, 1, 1, sin(state_A.yaw) * sin(state_A.pitch) * sin(state_A.roll) + cos(state_A.yaw) * cos(state_A.roll));
	gsl_matrix_set(M_TA, 1, 2, sin(state_A.yaw) * sin(state_A.pitch) * cos(state_A.roll) - cos(state_A.yaw) * sin(state_A.roll));
	gsl_matrix_set(M_TA, 1, 3, state_A.y);
	gsl_matrix_set(M_TA, 2, 0, -sin(state_A.pitch));
	gsl_matrix_set(M_TA, 2, 1, cos(state_A.pitch) * sin(state_A.roll));
	gsl_matrix_set(M_TA, 2, 2, cos(state_A.pitch) * cos(state_A.roll));
	gsl_matrix_set(M_TA, 2, 3, state_A.z);
	gsl_matrix_set(M_TA, 3, 0, 0.0);
	gsl_matrix_set(M_TA, 3, 1, 0.0);
	gsl_matrix_set(M_TA, 3, 2, 0.0);
	gsl_matrix_set(M_TA, 3, 3, 1.0);

	int signum;
	gsl_permutation * p = gsl_permutation_alloc(4);
	gsl_matrix* M_AT = gsl_matrix_alloc(4, 4);
	gsl_matrix* M_TA_temp = gsl_matrix_alloc(4, 4);
	gsl_matrix_memcpy(M_TA_temp, M_TA);
	gsl_linalg_LU_decomp (M_TA_temp, p, &signum);
	gsl_linalg_LU_invert (M_TA_temp, p, M_AT);

	// M_AB, transformation matrix for (state_A) -> (state_B)
	// M_AB = M_TB * (M_TA)^(-1)
	// M_AB = M_TB * M_AT
	gsl_matrix* M_AB = gsl_matrix_alloc(4, 4);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, M_TB, M_AT,  0.0, M_AB);

	gsl_permutation_free(p);
	gsl_matrix_free(M_TA);
	gsl_matrix_free(M_TA_temp);
	gsl_matrix_free(M_TB);
	gsl_matrix_free(M_AT);

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
	const double delta_Phi = state_B.roll - state_A.roll;
	const double delta_Theta = state_B.pitch - state_A.pitch;
	const double delta_Psi = state_B.yaw - state_A.yaw;

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

		const NeedleCoordinates intermediate_state = {
			state_A.x + intermediate_delta_X,
			state_A.y + intermediate_delta_Y,
			state_A.z + intermediate_delta_Z,
			state_A.roll + intermediate_delta_Phi,
			state_A.pitch + intermediate_delta_Theta,
			state_A.yaw + intermediate_delta_Psi
		};
		
		transformation_frames[j] = get_transformation_matrix(state_A, intermediate_state); 
	}
}

// Get vectors for sliders at their highest possible position (top_state).
Sliders_Model_Nullable* get_sliders_in_top_state(const double radius, const double angle_in_radians, bool print_debug_info_local) {
	Sliders_Model_Nullable* sliders = (Sliders_Model_Nullable*) malloc(sizeof(Sliders_Model_Nullable));

	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
		sliders->vecs[i] = gsl_vector_alloc(3);
	}

	gsl_vector_set(sliders->vecs[0], 0, radius * sin(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->vecs[0], 1, radius * cos(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->vecs[0], 2, MAX_SLIDER_HEIGHT_CM);
	gsl_vector_set(sliders->vecs[1], 0, radius * sin(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->vecs[1], 1, radius * cos(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->vecs[1], 2, MAX_SLIDER_HEIGHT_CM);
	gsl_vector_set(sliders->vecs[2], 0, radius * sin(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->vecs[2], 1, radius * cos(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->vecs[2], 2, MAX_SLIDER_HEIGHT_CM);
	gsl_vector_set(sliders->vecs[3], 0, radius * sin(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->vecs[3], 1, radius * cos(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->vecs[3], 2, MAX_SLIDER_HEIGHT_CM);
	gsl_vector_set(sliders->vecs[4], 0, radius * sin(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->vecs[4], 1, radius * cos(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->vecs[4], 2, MAX_SLIDER_HEIGHT_CM);
	gsl_vector_set(sliders->vecs[5], 0, radius * sin(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->vecs[5], 1, radius * cos(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(sliders->vecs[5], 2, MAX_SLIDER_HEIGHT_CM);


	if (print_debug_info_local) {
		printf("\n==== Slider coordinates (x, y, z) in cm for (top_state):\n");
		print_sliders(sliders);
	}

	return sliders;
}

Table_Model* get_table_in_top_state(bool print_debug_info_local) {
	Table_Model* table = (Table_Model*)malloc(sizeof(Table_Model));
	for (size_t i = 0; i < NUM_OF_VECTORS; i++) {
		table->vecs[i] = gsl_vector_alloc(3);
	}

	Sliders_Model_Nullable* sliders_in_top_state = get_sliders_in_top_state(SLIDERS_RADIUS_CM, SLIDERS_ANGLE_RAD, false);

	double radius = TABLE_RADIUS_CM;
	double angle_in_radians = TABLE_ANGLE_RAD;
	gsl_vector_set(table->vecs[0], 0, radius * sin(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(table->vecs[0], 1, radius * cos(DEGREES_TO_RADIANS * 0 + (angle_in_radians / 2.0)));
	gsl_vector_set(table->vecs[1], 0, radius * sin(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(table->vecs[1], 1, radius * cos(DEGREES_TO_RADIANS * 0 - (angle_in_radians / 2.0)));
	gsl_vector_set(table->vecs[2], 0, radius * sin(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(table->vecs[2], 1, radius * cos(DEGREES_TO_RADIANS * 120 + (angle_in_radians / 2.0)));
	gsl_vector_set(table->vecs[3], 0, radius * sin(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(table->vecs[3], 1, radius * cos(DEGREES_TO_RADIANS * 120 - (angle_in_radians / 2.0)));
	gsl_vector_set(table->vecs[4], 0, radius * sin(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(table->vecs[4], 1, radius * cos(DEGREES_TO_RADIANS * 240 + (angle_in_radians / 2.0)));
	gsl_vector_set(table->vecs[5], 0, radius * sin(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));
	gsl_vector_set(table->vecs[5], 1, radius * cos(DEGREES_TO_RADIANS * 240 - (angle_in_radians / 2.0)));

	// Now we calculate the z-coordinate for the first table vertex. For (top_state) all 6 hexagon vertices and centroid will have same z-coordinate.
	double slider_x = gsl_vector_get(sliders_in_top_state->vecs[0], 0);
	double slider_y = gsl_vector_get(sliders_in_top_state->vecs[0], 1);
	double slider_z = gsl_vector_get(sliders_in_top_state->vecs[0], 2);
	double table_x = gsl_vector_get(table->vecs[0], 0);
	double table_y = gsl_vector_get(table->vecs[0], 1);
	double deltaX = slider_x - table_x;
	double deltaY = slider_y - table_y;
	double deltaZ_squared = pow(SLIDER_ARM_LENGTH_CM, 2) - pow(deltaX, 2) - pow(deltaY, 2);
	assert(deltaZ_squared > 0);
	double table_z = slider_z - sqrt(deltaZ_squared); 

	for (size_t i = 0; i < NUM_OF_HEXAGON_VERTICES; i++) {
		gsl_vector_set(table->vecs[i], 2, table_z);
	}

	gsl_vector_set(table->vecs[6], 0, 0.0); 
	gsl_vector_set(table->vecs[6], 1, 0.0);
	gsl_vector_set(table->vecs[6], 2, table_z); 

	gsl_vector_set(table->vecs[7], 0, 0.0);
	gsl_vector_set(table->vecs[7], 1, 0.0);
	gsl_vector_set(table->vecs[7], 2, table_z - NEEDLE_LENGTH_CM);

	if (print_debug_info_local) {
		printf("\n==== Table coordinates (x, y, z) in cm for (top_state):\n");
		print_table(table);
	}

	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
	    gsl_vector_free(sliders_in_top_state->vecs[i]);
	}
	free(sliders_in_top_state);

	return table;
}

Sliders_Model_Nullable* get_sliders(const Table_Model table_transformed_vecs, const Sliders_Model_Nullable sliders_in_top_state) {
	Sliders_Model_Nullable* sliders = (Sliders_Model_Nullable*) malloc(sizeof(Sliders_Model_Nullable));

	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
		double slider_x = gsl_vector_get(sliders_in_top_state.vecs[i], 0);
		double slider_y = gsl_vector_get(sliders_in_top_state.vecs[i], 1);

		double deltaX = gsl_vector_get(sliders_in_top_state.vecs[i], 0) - gsl_vector_get(table_transformed_vecs.vecs[i], 0);
		double deltaY = gsl_vector_get(sliders_in_top_state.vecs[i], 1) - gsl_vector_get(table_transformed_vecs.vecs[i], 1);
		double deltaZ_squared = pow(SLIDER_ARM_LENGTH_CM, 2) - pow(deltaX, 2) - pow(deltaY, 2);

		if (deltaZ_squared > 0) {
			// This is one of two possible solutions for z coordinate of slider, the higher one;
			// Lower solution will have minus sqrt(deltaZ_squared) instead of plus sqrt(deltaZ_squared);
			// We don't use second solution at all because we fear damage during possible solution switching.
			double slider_z = gsl_vector_get(table_transformed_vecs.vecs[i], 2) + sqrt(deltaZ_squared);
			if ((slider_z <= MAX_SLIDER_HEIGHT_CM) && (slider_z >= MIN_SLIDER_HEIGHT_CM)) {
				sliders->vecs[i] = gsl_vector_alloc(3);
				gsl_vector_set(sliders->vecs[i], 0, slider_x);
				gsl_vector_set(sliders->vecs[i], 1, slider_y);
				gsl_vector_set(sliders->vecs[i], 2, slider_z);
			} else {
				sliders->vecs[i] = NULL;
			};
		} else {
			sliders->vecs[i] = NULL;
		}
	}
	return sliders;
}

int set_sliders_frames(
		Sliders_Model_Nullable* slider_frames,
		const size_t num_of_frames,
		const Table_Model* table_frames,
		const Sliders_Model_Nullable sliders_in_top_state
		) {

	for (size_t j = 0; j < num_of_frames; j++)  {
		Table_Model table_for_this_frame = table_frames[j];
		Sliders_Model_Nullable* sliders_for_this_frame = get_sliders(table_for_this_frame, sliders_in_top_state);
		
		for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
			if (sliders_for_this_frame->vecs[i] == 0) {
				printf("\n======== FAILED to compute sliders\n");
				return 1;
			}
		}

		slider_frames[j] = *sliders_for_this_frame;
	}
	return 0;
}

Table_Model* get_transformed_table(const gsl_matrix* T, const Table_Model initial_table){
	Table_Model* transformed_table = (Table_Model*)malloc(sizeof(Table_Model));

	for (size_t i = 0; i < NUM_OF_VECTORS; i++) {
		const gsl_vector* initial_vec_in_3D = initial_table.vecs[i];

		// Uplifting to 4D
		gsl_vector* initial_vec_in_4D = gsl_vector_alloc(4);
		gsl_vector_set(initial_vec_in_4D, 0, gsl_vector_get(initial_vec_in_3D, 0)); // We take X coordinate of the 3D vector
		gsl_vector_set(initial_vec_in_4D, 1, gsl_vector_get(initial_vec_in_3D, 1)); // We take Y coordinate of the 3D vector
		gsl_vector_set(initial_vec_in_4D, 2, gsl_vector_get(initial_vec_in_3D, 2)); // We take Z coordinate of the 3D vector
		gsl_vector_set(initial_vec_in_4D, 3, 1.0);                                  // We use 1.0 as coordinate for 4th dim

		// Linear transformation in 4D
		gsl_vector* transformed_vec_in_4D = gsl_vector_alloc(4);
		gsl_blas_dgemv(CblasNoTrans, 1.0, T, initial_vec_in_4D, 0.0, transformed_vec_in_4D);

		/*// Projection back to 3D*/
		transformed_table->vecs[i] = gsl_vector_alloc(3);
		gsl_vector_set(transformed_table->vecs[i], 0, gsl_vector_get(transformed_vec_in_4D, 0));
		gsl_vector_set(transformed_table->vecs[i], 1, gsl_vector_get(transformed_vec_in_4D, 1));
		gsl_vector_set(transformed_table->vecs[i], 2, gsl_vector_get(transformed_vec_in_4D, 2));

		gsl_vector_free(initial_vec_in_4D);
		gsl_vector_free(transformed_vec_in_4D);
	}
	return transformed_table;
}

NeedleCoordinates* get_top_state(bool print_debug_info_local) {
	Table_Model* table_in_top_state = get_table_in_top_state(false);
	double needle_z_in_top_state = gsl_vector_get(table_in_top_state->vecs[7], 2);

	NeedleCoordinates* top_state = (NeedleCoordinates*) malloc(sizeof(NeedleCoordinates));
	*top_state = (NeedleCoordinates) {0.0, 0.0, needle_z_in_top_state, 0.0, 0.0, 0.0};

	if (print_debug_info_local) {
		printf("==== Needle coordinates for (top_state):\n");
		print_needle_coordinates(top_state);
	}

	for (size_t i = 0; i < NUM_OF_VECTORS; i++) {
	    gsl_vector_free(table_in_top_state->vecs[i]);
	}
	free(table_in_top_state);

	return top_state;
}

Trajectory* get_trajectory(const NeedleCoordinates initial_state, const NeedleCoordinates final_state, uint64_t movement_start_mks, const size_t num_of_frames) {
	const uint64_t movement_duration_mks = MOVEMENT_DURATION_MKS; // TODO? - calculate this
	assert(num_of_frames >= 2);

	NeedleCoordinates* top_state = get_top_state(PRINT_DEBUG_INFO_GLOBAL);
	Table_Model* table_in_top_state = get_table_in_top_state(PRINT_DEBUG_INFO_GLOBAL);
	gsl_matrix* M_TA = get_transformation_matrix(*top_state, initial_state);
	Table_Model* table_in_state_A = get_transformed_table(M_TA, *table_in_top_state);
	Sliders_Model_Nullable* sliders_in_top_state = get_sliders_in_top_state(SLIDERS_RADIUS_CM, SLIDERS_ANGLE_RAD, PRINT_DEBUG_INFO_GLOBAL);

	// First frame is the the identity transformation (Id).
	// Last frame is the full transformation (M_AB), which brings (state_A) to the (state_B).
	// Other frames describe transformation that is "in progress", which bring (state_A) to some intermediate states
	// Example for the case with 4 frames:
	// 	frame 0 is the Id;
	//	frame 1 is some mix between Id and M_AB (less M_AB)
	//	frame 2 is some mix between Id and M_AB (more M_AB)
	//	frame 3 is the M_AB;
	gsl_matrix* M_frames[num_of_frames];
	set_transformation_frames(M_frames, initial_state, final_state, num_of_frames);

	// Each 3D model is produced by applying one of transformations "in progress" to the 3D model for initial state.
	Table_Model table_frames[num_of_frames];
	for (size_t j = 0; j < num_of_frames; j++) {
		table_frames[j] = *get_transformed_table(M_frames[j], *table_in_state_A);
	}

	Sliders_Model_Nullable slider_frames[num_of_frames];
	int status_code = set_sliders_frames(slider_frames, num_of_frames, table_frames, *sliders_in_top_state);
	if (status_code == 1) {
		printf("\n======== FAILED to create trajectory between (initial_state) and (final_state). Please check if states are set correctly:");
		printf("\n==== Needle coordinates for the (initial_state):\n");
		print_needle_coordinates(&initial_state);
		printf("\n==== Needle coordinates for the (final_state):\n");
		print_needle_coordinates(&final_state);
		return NULL;
	}

	uint64_t timestamps[num_of_frames];
	for (size_t i = 0; i < num_of_frames; i++) {
		const double progress = (double) i / (double)(num_of_frames - 1);
	    	timestamps[i] = movement_start_mks + (uint64_t)(movement_duration_mks * progress);
	}

	Trajectory* trajectory = (Trajectory*)malloc(sizeof(Trajectory));
	trajectory->initial_state = malloc(sizeof(NeedleCoordinates));
	trajectory->final_state = malloc(sizeof(NeedleCoordinates));

	for (size_t i = 0; i < num_of_frames; i++) {
		trajectory->table_frames[i] = table_frames[i];
		trajectory->sliders_frames[i] = slider_frames[i];
		trajectory->timestamps[i] = timestamps[i];
	}

	memcpy(trajectory->initial_state, &initial_state, sizeof(NeedleCoordinates));
	memcpy(trajectory->final_state, &final_state, sizeof(NeedleCoordinates));
	trajectory->num_of_frames = num_of_frames;

	if (PRINT_DEBUG_INFO_GLOBAL) {
	    print_trajectory(trajectory, true);
	}

	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
	    gsl_vector_free(sliders_in_top_state->vecs[i]);
	}
	free(sliders_in_top_state);
	for (size_t i = 0; i < NUM_OF_VECTORS; i++) {
	    gsl_vector_free(table_in_top_state->vecs[i]);
	}
	free(table_in_top_state);
	for (size_t i = 0; i < NUM_OF_VECTORS; i++) {
	    gsl_vector_free(table_in_state_A->vecs[i]);
	}
	free(table_in_state_A);
	for (size_t j = 0; j < num_of_frames; j++) {
	        gsl_matrix_free(M_frames[j]);
	}
	gsl_matrix_free(M_TA);
	free(top_state);

	return trajectory;
}

uint64_t get_current_time_in_mks() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000 * 1000 + (uint64_t)ts.tv_nsec / 1000;
}

void free_trajectory(Trajectory* trajectory) {
    for (size_t i = 0; i < trajectory->num_of_frames; i++) {
        for (size_t j = 0; j < NUM_OF_VECTORS; j++) {
            gsl_vector_free(trajectory->table_frames[i].vecs[j]);
        }
    }
    for (size_t i = 0; i < trajectory->num_of_frames; i++) {
        for (size_t j = 0; j < NUM_OF_SLIDERS; j++) {
            if (trajectory->sliders_frames[i].vecs[j] != NULL) {
                gsl_vector_free(trajectory->sliders_frames[i].vecs[j]);
            }
        }
    }
    free(trajectory->initial_state);
    free(trajectory->final_state);
    free(trajectory);
}

// Python - C interface
bool set_frames_py(
		double* table_x,
		double* table_y,
		double* table_z,
		double* slider_x,
		double* slider_y,
		double* slider_z,
		const double* ptr_state_A,
		const double* ptr_state_B,
		const size_t requested_num_of_frames
){
	const NeedleCoordinates state_A = {ptr_state_A[0], ptr_state_A[1], ptr_state_A[2], ptr_state_A[3], ptr_state_A[4], ptr_state_A[5]};
	const NeedleCoordinates state_B = {ptr_state_B[0], ptr_state_B[1], ptr_state_B[2], ptr_state_B[3], ptr_state_B[4], ptr_state_B[5]};

	uint64_t movement_start_mks = get_current_time_in_mks() + 100;
	Trajectory* trajectory = get_trajectory(state_A, state_B, movement_start_mks, requested_num_of_frames);
    
	if (trajectory == NULL) {
		return false; 
	}
	bool success = true;
	for (size_t i = 0; i < requested_num_of_frames; ++i) {
	    for (int j = 0; j < NUM_OF_VECTORS; ++j) {
		table_x[i * NUM_OF_VECTORS + j] = gsl_vector_get(trajectory->table_frames[i].vecs[j], 0);
		table_y[i * NUM_OF_VECTORS + j] = gsl_vector_get(trajectory->table_frames[i].vecs[j], 1);
		table_z[i * NUM_OF_VECTORS + j] = gsl_vector_get(trajectory->table_frames[i].vecs[j], 2);
	    }
	    for (int k = 0; k < NUM_OF_SLIDERS; ++k) {
		if (trajectory->sliders_frames[i].vecs[k] != NULL) {
		    slider_x[i * NUM_OF_SLIDERS + k] = gsl_vector_get(trajectory->sliders_frames[i].vecs[k], 0);
		    slider_y[i * NUM_OF_SLIDERS + k] = gsl_vector_get(trajectory->sliders_frames[i].vecs[k], 1);
		    slider_z[i * NUM_OF_SLIDERS + k] = gsl_vector_get(trajectory->sliders_frames[i].vecs[k], 2);
		} else {
		    slider_x[i * NUM_OF_SLIDERS + k] = NAN;
		    slider_y[i * NUM_OF_SLIDERS + k] = NAN;
		    slider_z[i * NUM_OF_SLIDERS + k] = NAN;
		    success = false;
		}
	    }
	}
	free_trajectory(trajectory);

	return success;
}


size_t get_next_timestamp_index(const uint64_t input_timestamp, const size_t num_of_timestamps, const uint64_t* timestamps, const bool print_debug_info_local) {
    size_t result_index = SIZE_MAX; // SIZE_MAX indicates failure in timestamp search

    if (input_timestamp < timestamps[0]) {
	result_index = 0;
    }
    else if (input_timestamp > timestamps[num_of_timestamps - 1]) {
        result_index = SIZE_MAX;
    } else {
        // Binary search to find exact or next timestamp
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
            printf("\nSearching for timestamp %llu microseconds\nNext timestamp found:   %llu microseconds (frame %lu / %lu)\n",
			    input_timestamp, timestamps[result_index], result_index, num_of_timestamps - 1);
        } else {
            printf("\nSearching for timestamp %llu microseconds\nFAILED to find next timestamp found\n", input_timestamp);
        }
    }

    return result_index;
}

Sliders_Movement_Recipe compare_with_trajectory(const Trajectory* trajectory, const Sliders_Model_Nullable* sliders_inferred, const size_t frame_index, bool print_debug_info_local){
	Sliders_Model_Nullable sliders_at_trajectory = trajectory->sliders_frames[frame_index];
	Sliders_Model_Nullable difference = get_slider_difference(sliders_at_trajectory, *sliders_inferred);
	Sliders_Movement_Recipe sliders_movement_recipe = {
		.slider_movement_needed[0] = fabs(gsl_vector_get(difference.vecs[0], 2)) >= MIN_SLIDER_STEP_CM, 
		.slider_movement_needed[1] = fabs(gsl_vector_get(difference.vecs[1], 2)) >= MIN_SLIDER_STEP_CM, 
		.slider_movement_needed[2] = fabs(gsl_vector_get(difference.vecs[2], 2)) >= MIN_SLIDER_STEP_CM, 
		.slider_movement_needed[3] = fabs(gsl_vector_get(difference.vecs[3], 2)) >= MIN_SLIDER_STEP_CM, 
		.slider_movement_needed[4] = fabs(gsl_vector_get(difference.vecs[4], 2)) >= MIN_SLIDER_STEP_CM, 
		.slider_movement_needed[5] = fabs(gsl_vector_get(difference.vecs[5], 2)) >= MIN_SLIDER_STEP_CM, 
		.slider_movement_direction_is_up[0] = gsl_vector_get(difference.vecs[0], 2) >= 0.0,
		.slider_movement_direction_is_up[1] = gsl_vector_get(difference.vecs[1], 2) >= 0.0,
		.slider_movement_direction_is_up[2] = gsl_vector_get(difference.vecs[2], 2) >= 0.0,
		.slider_movement_direction_is_up[3] = gsl_vector_get(difference.vecs[3], 2) >= 0.0,
		.slider_movement_direction_is_up[4] = gsl_vector_get(difference.vecs[4], 2) >= 0.0,
		.slider_movement_direction_is_up[5] = gsl_vector_get(difference.vecs[5], 2) >= 0.0,
	};

	if (PRINT_DEBUG_INFO_GLOBAL && print_debug_info_local) {
		printf("\n=== Comparing current sliders and sliders expected for frame %lu / %lu", frame_index, trajectory->num_of_frames -  1);
		printf("\nCurrent sliders (inferred from past inputs to our physical system):\n");
		print_sliders(sliders_inferred);
		printf("\nExpected sliders:\n");
		print_sliders(&sliders_at_trajectory);
		printf("\nDifference in sliders:\n");
		print_sliders(&difference);
		print_sliders_movement_recipe(&sliders_movement_recipe);
	}

	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
	    if (difference.vecs[i] != NULL) {
		gsl_vector_free(difference.vecs[i]);
	    }
	}

	return sliders_movement_recipe;
}


void mutate_infered_slider_state (Sliders_Model_Nullable* current_sliders, const Sliders_Movement_Recipe sliders_movement_recipe) {
	for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
		if (sliders_movement_recipe.slider_movement_needed[i]) {
			if (sliders_movement_recipe.slider_movement_direction_is_up[i]) { 
				gsl_vector_set(current_sliders->vecs[i], 2, gsl_vector_get(current_sliders->vecs[i], 2) + MIN_SLIDER_STEP_CM);
			} else {
				gsl_vector_set(current_sliders->vecs[i], 2, gsl_vector_get(current_sliders->vecs[i], 2) - MIN_SLIDER_STEP_CM);
			}
		} 
	}
}

bool trajectory_converged (Sliders_Movement_Recipe recipe, size_t index_of_frame) {
	return 				
		!recipe.slider_movement_needed[0] &&
		!recipe.slider_movement_needed[1] &&
		!recipe.slider_movement_needed[2] &&
		!recipe.slider_movement_needed[3] &&
		!recipe.slider_movement_needed[4] &&
		index_of_frame == NUM_OF_FRAMES - 1;
}

int follow_trajectory(Sliders_Model_Nullable* current_sliders, const Trajectory* trajectory) {
	const size_t num_of_frames = trajectory->num_of_frames;
	size_t num_of_iterations = 50 * 1000;
	size_t sleep_duration_mks = MOVEMENT_DURATION_MKS / num_of_iterations;
	size_t debug_print_interval = (num_of_iterations + NUM_OF_FRAMES_TO_PRINT_ON_FOLLOWING - 1) / NUM_OF_FRAMES_TO_PRINT_ON_FOLLOWING; // This is for debug info printing only.

	if (PRINT_DEBUG_INFO_GLOBAL) {
		printf("\n======== Starting this trajectory:");
		print_trajectory(trajectory, false);
	}

	for (size_t i = 0; i <= num_of_iterations; i++) {
		bool print_debug_info_local = (i % debug_print_interval == 0);
		size_t index_of_the_next_frame = get_next_timestamp_index(get_current_time_in_mks(), num_of_frames, trajectory->timestamps, print_debug_info_local);

		// If (index_of_the_next_frame != SIZE_MAX) then next frame was found
		if (index_of_the_next_frame != SIZE_MAX) {
			Sliders_Movement_Recipe recipe = compare_with_trajectory(trajectory, current_sliders, index_of_the_next_frame, print_debug_info_local);
			mutate_infered_slider_state(current_sliders, recipe);
			if (trajectory_converged(recipe, index_of_the_next_frame)) {

				if (PRINT_DEBUG_INFO_GLOBAL) {
					printf("\n======== This trajectory was completed sucessfully:");
					print_trajectory(trajectory, false);
				}
				return 0;
			}
		}
		usleep(sleep_duration_mks);
	} 
	printf("\n======== FAILED to complete trajectory");
	return 1;
}

Sliders_Model_Nullable* deepcopy_sliders(const Sliders_Model_Nullable* original) {
    Sliders_Model_Nullable* deepcopy = malloc(sizeof(Sliders_Model_Nullable));
    for (size_t i = 0; i < NUM_OF_SLIDERS; i++) {
        if (original->vecs[i] != NULL) {
            deepcopy->vecs[i] = gsl_vector_alloc(original->vecs[i]->size);
            gsl_vector_memcpy(deepcopy->vecs[i], original->vecs[i]);
        } else {
            deepcopy->vecs[i] = NULL;
        }
    }
    return deepcopy;
}

int main(void) {
	// Hardware interaction simulation

	// Moving to the (top_state). This is the state with all sliders at top.
	const NeedleCoordinates initial_state = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};	// Assumed initial state
	const size_t num_of_frames = NUM_OF_FRAMES; // TODO? - calculate this
	const NeedleCoordinates top_state = *get_top_state(false);
	Trajectory* trajectory_0 = get_trajectory(initial_state, top_state, get_current_time_in_mks() + 100, num_of_frames);
	assert(trajectory_0 != NULL); // Check if trajectory was created successfully
	Sliders_Model_Nullable* sliders = deepcopy_sliders(&trajectory_0->sliders_frames[0]); // Sliders represent history of inputs to the physical system, do not throw them away. 
	const int status_code_0 = follow_trajectory(sliders, trajectory_0);
	assert(status_code_0 == 0); // Check if trajectory was followed successfully
	free_trajectory(trajectory_0);

	// Moving from the (top_state) to the (A_state):
	const NeedleCoordinates A_state = {-2.58, 2.77, 3.94, 0.71, 1.22, 0.47};
	Trajectory* trajectory_1 = get_trajectory(top_state, A_state, get_current_time_in_mks() + 100, num_of_frames);
	assert(trajectory_1 != NULL);
	const int status_code_1 = follow_trajectory(sliders, trajectory_1);
	assert(status_code_1 == 0);
	free_trajectory(trajectory_1);

	// Moving back to the (top_state):
	Trajectory* trajectory_2 = get_trajectory(A_state, top_state, get_current_time_in_mks() + 100, num_of_frames);
	assert(trajectory_2 != NULL);
	const int status_code_2 = follow_trajectory(sliders, trajectory_2);
	assert(status_code_2 == 0);
	free_trajectory(trajectory_2);

	printf("\n======== Success");
	return 0;
}

