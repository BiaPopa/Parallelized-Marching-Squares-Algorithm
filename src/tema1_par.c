// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }
#define MIN(a, b) ((a) < (b) ? (a) : (b))

struct thread_args {
    unsigned char **grid;
    long id;
    long no_of_threads;
    int step_x;
    int step_y;
    int start;	
    int end;
    ppm_image *image;
    ppm_image *image_aux;
    ppm_image *rescaled_image;
    ppm_image **contour_map;
    pthread_barrier_t *barrier;
};

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(struct thread_args *argument, unsigned char k, int x_i, int y_j) {	
    ppm_image *contour = argument->contour_map[k];
    int x = x_i * argument->step_x;
    int y = y_j * argument->step_y;

    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * argument->rescaled_image->y + y + j;

            argument->rescaled_image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            argument->rescaled_image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            argument->rescaled_image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
void sample_grid_thread(struct thread_args *argument, unsigned char sigma) {
    int p = argument->rescaled_image->x / argument->step_x;
    int q = argument->rescaled_image->y / argument->step_y;
    
    // Only the first thread allocates memory for the grid
    // The other threads will point to the same memory address
    if (argument->id == 0) {
        argument->grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
        if (!argument->grid) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }

        for (int i = 0; i <= p; i++) {
            argument->grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
            if (!argument->grid[i]) {
                fprintf(stderr, "Unable to allocate memory\n");
                exit(1);
            }
        }

        for (long id = 0; id < argument->no_of_threads; id++) {
            argument[id].grid = argument->grid;
        }
    }

    pthread_barrier_wait(argument->barrier);

    // Compute the start and end indices using the formula from lab 1
    argument->start = argument->id * p / argument->no_of_threads;
    argument->end = MIN((argument->id + 1) * p / argument->no_of_threads, p);

    // Each thread iterates over a part of the image
    for (int i = argument->start; i < argument->end; i++) {	
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = argument->rescaled_image->data[i * argument->step_x * argument->rescaled_image->y + j * argument->step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > sigma) {
                argument->grid[i][j] = 0;
            } else {
                argument->grid[i][j] = 1;
            }
        }
    }
    argument->grid[p][q] = 0;

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    // Each thread iterates over a part of the image
    for (int i = argument->start; i < argument->end; i++) {
        ppm_pixel curr_pixel = argument->rescaled_image->data[i * argument->step_x * argument->rescaled_image->y + argument->rescaled_image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            argument->grid[i][q] = 0;
        } else {
            argument->grid[i][q] = 1;
        }
    }

    // Compute the start and end indices using the formula from lab 1
    argument->start = argument->id * q / argument->no_of_threads;
    argument->end = MIN((argument->id + 1) * q / argument->no_of_threads, q);

    // Each thread iterates over a part of the image
    for (int j = argument->start; j < argument->end; j++) {
        ppm_pixel curr_pixel = argument->rescaled_image->data[(argument->rescaled_image->x - 1) * argument->rescaled_image->y + j * argument->step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            argument->grid[p][j] = 0;
        } else {
            argument->grid[p][j] = 1;
        }
    }
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march_thread(struct thread_args *argument) {
    int p = argument->rescaled_image->x / argument->step_x;
    int q = argument->rescaled_image->y / argument->step_y;

    // Compute the start and end indices using the formula from lab 1
    argument->start = argument->id * p / argument->no_of_threads;
    argument->end = MIN((argument->id + 1) * p / argument->no_of_threads, p);

    // Each thread iterates over a part of the image
    for (int i = argument->start; i < argument->end; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * argument->grid[i][j] + 4 * argument->grid[i][j + 1] + 2 * argument->grid[i + 1][j + 1] + 1 * argument->grid[i + 1][j];
            update_image(argument, k, i, j);
        }
    }
}

void rescale_image_thread(struct thread_args *argument) {
    uint8_t sample[3];

    // Only the first thread allocates memory for the auxiliary image
    // The other threads will point to the same memory address
    if (argument->id == 0) {
        argument->image_aux = (ppm_image*)malloc(sizeof(ppm_image));
        if (!argument->image_aux) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
        argument->image_aux->x = RESCALE_X;
        argument->image_aux->y = RESCALE_Y;

        argument->image_aux->data = (ppm_pixel*)malloc(argument->image_aux->x * argument->image_aux->y * sizeof(ppm_pixel));
        if (!argument->image_aux) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }

        for (long id = 1; id < argument->no_of_threads; id++) {
            argument[id].image_aux = argument->image_aux;
        }
    }

    pthread_barrier_wait(argument->barrier);

    // Compute the start and end indices using the formula from lab 1
    argument->start = argument->id * argument->image_aux->x / argument->no_of_threads;
    argument->end = MIN((argument->id + 1) * argument->image_aux->x / argument->no_of_threads, argument->image_aux->x);

    // Each thread iterates over a part of the image
    for (int i = argument->start; i < argument->end; i++) {
        for (int j = 0; j < argument->image_aux->y; j++) {
            float u = (float)i / (float)(argument->image_aux->x - 1);
            float v = (float)j / (float)(argument->image_aux->y - 1);
            sample_bicubic(argument->image, u, v, sample);

            argument->image_aux->data[i * argument->image_aux->y + j].red = sample[0];
            argument->image_aux->data[i * argument->image_aux->y + j].green = sample[1];
            argument->image_aux->data[i * argument->image_aux->y + j].blue = sample[2];
        }
    }

    pthread_barrier_wait(argument->barrier);

    // Only the first thread frees the memory of the initial image
    if (argument->id == 0) {
        free(argument->image->data);
        free(argument->image);
    }
    
    // rescaled_image is image_aux after all threads have finished their rescaling
    argument->rescaled_image = argument->image_aux;
}

void *thread_func(void *arg) {
    struct thread_args *arguments = (struct thread_args *)arg;
    
    // 1. Rescale the image
    if (arguments->image->x <= RESCALE_X && arguments->image->y <= RESCALE_Y) {
        arguments->rescaled_image = arguments->image;
    } else {
        rescale_image_thread(arguments);
        pthread_barrier_wait(arguments->barrier);
    }

    // 2. Sample the grid
    sample_grid_thread(arguments, SIGMA);
    pthread_barrier_wait(arguments->barrier);

    // 3. March the squares
    march_thread(arguments);
    pthread_barrier_wait(arguments->barrier);

    pthread_exit(NULL);
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    long no_of_threads = atoi(argv[3]);
    pthread_t threads[no_of_threads];
    struct thread_args arguments[no_of_threads];
    ppm_image *scaled_image;

    unsigned char **grid;
    int step_x = STEP;
    int step_y = STEP;
    int r;

    // Initialize barrier
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, no_of_threads);

    // Read image
    ppm_image *image = read_ppm(argv[1]);

    // 0. Initialize contour map
    ppm_image **contour_map = init_contour_map();

    // Create no_of_threads threads
    for (long id = 0; id < no_of_threads; id++) {
        arguments[id].barrier = &barrier;
        arguments[id].id = id;
        arguments[id].no_of_threads = no_of_threads;
        arguments[id].image = image;
        arguments[id].step_x = step_x;
        arguments[id].step_y = step_y;
        arguments[id].contour_map = contour_map;

        r = pthread_create(&threads[id], NULL, thread_func, &arguments[id]);

        if (r) {
            printf("Eroare la crearea thread-ului %ld\n", id);
            exit(-1);
        }
    }

    // Wait for all threads to finish their execution
    for (long id = 0; id < no_of_threads; id++) {
        r = pthread_join(threads[id], NULL);

        if (r) {
            printf("Eroare la asteptarea thread-ului %ld\n", id);
            exit(-1);
        }
    }

    // Save the final results from the structure of the first thread
    scaled_image = arguments[0].rescaled_image;
    grid = arguments[0].grid;

    // 4. Write output
    write_ppm(scaled_image, argv[2]);

    // Destroy barrier and free all resources
    pthread_barrier_destroy(&barrier);
    free_resources(scaled_image, contour_map, grid, step_x);
    return 0;
}
