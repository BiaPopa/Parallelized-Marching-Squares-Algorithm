Popa Bianca

Parallelized Marching Squares Algorithm

In main, I extracted the number of threads from argv[3] and declared
an array of pthread_t and one of struct thread_args. The structure
contains the following fields:

unsigned char **grid;        -- grid of the given image \
long id;                     -- thread id \
long no_of_threads;          -- number of threads \
int step_x;                  -- step for x (8) \
int step_y;                  -- step for y (8)        
int start;	                 -- start index \
int end;                     -- end index \
ppm_image *image;            -- initial image \
ppm_image *image_aux;        -- auxiliary image (used for rescaling) \
ppm_image *rescaled_image;   -- rescaled image \
ppm_image **contour_map;     -- contour map (from init_contour_map) \
pthread_barrier_t *barrier;  -- barrier 

I initialized the barrier and created no_of_threads threads. For each thread, 
I initialized their structure and called the pthread_create function. Then, 
waited for all threads to finish their execution using pthread_join.
Finally, I destroyed the barrier and freed all the allocated memory.

I moved the functionality of Marching Squares algorithm in the thread_func:

1. Rescale the image

If the image has a smaller size than RESCALE_X and RESCALE_Y, then rescaled_image
will be the initial image. Otherwise, rescale_image_thread function is called.
It has the same functionality as rescale_image(ppm_image *image) from secvential
solution, but with the following differences:

Only the first thread will alocate memory for the auxiliary image, the others
will point to the same memory address. Before rescaling the initial image,
the program will compute the start and end indices using the formula:

int start = ID * N / P; \
int end = min((ID + 1) * N / P, N);

where N is the number of elements we want to iterate through and P is the number
of threads. Then, each thread will iterate through their own interval (from start
to end).

The program waits for all threads to finish rescaling the image. The first thread
will free the memory allocated for the initial image and all threads will initialize
rescaled_image with the auxiliary image.

2. Sample the grid

The program calls sample_grid_thread function, which has the same functionality as
sample_grid(ppm_image *image, int step_x, int step_y, unsigned char sigma), but with
the following differences:

Only the first thread will alocate memory for the grid, the others will point to the
same memory address. The program will compute the start and end indices using the 
formula from above. Then, each thread will iterate through their own interval (from
start to end). The first two iterations are divided after p (argument->rescaled_image
->x / argument->step_x) and the last one after q (argument->rescaled_image->y / 
argument->step_y).

3. March the squares

Finally, thread_func calls march_thread function, which has the same functionality as
march(ppm_image *image, unsigned char **grid, ppm_image **contour_map, int step_x, int 
step_y), but with the following differences:

The program will compute the start and end indices using the formula from first step.
Then, each thread will iterate through their own interval (from start to end). 
update_image was also modified to work with the structure of the threads.

init_contour_map and free_resources functions remain unchanged.

Test 6 speedup: 2.29 for 2 threads, 4.17 for 4 threads \
Test 7 speedup: 2.42 for 2 threads, 4.70 for 4 threads

