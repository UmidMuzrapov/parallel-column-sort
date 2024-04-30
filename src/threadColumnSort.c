/**
 * @author Umid Muzrapov
 * Instructor: David Lowenthal
 *
 * Description: The program implements the parallelized version of a column sort function
 * that sorts nonnegative integers using pthreads. The full description of
 * the algorithm can be found in the requirements doc.
 *
 * Operational Requirements:
 *  C99
 *  stdio.h
 *  stdlib.h
 *  pthread.h
 *  math.h
 *  sys/time.h
 *  columSort.h
 *
 */

#include "columnSort.h"
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>
#include <semaphore.h>

/***** Macros and Type Definitions*****/
#define MAX 1000000
#define MIN (-1000000)

typedef struct
{
    int start;
    int end;
    int id;
} worker_data_t;

typedef struct
{
    volatile int *arrive;
    volatile int done;
} dissemination_barrier_t;

/***** Global Variables *****/
int size, NUMBER_OF_THREADS, LENGTH, work_size, number_of_stages;
int **matrix;
dissemination_barrier_t barrier;

// semaphores used to turn between sorting and nonsorting steps.
sem_t nonsorting_sem, sorting_sem;
// task complete singles that sorting threads can break the infinite loop
int task_complete = 0;

/***** Utility Functions *****/

/**
 * print_matrix function prints the matrix and
 * is used for debugging purposes.
 *
 * @param twoDArray int** used to point the rows of the array.
 * @param step debugging step
 * @param rows number of rows in 2-d matrix
 * @param cols number of columns in 2-d matrix
 */
void printMatrix(int step, int rows, int cols)
{
    printf("------------%d\n", step);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%10d ", matrix[i][j]);
        }

        printf("\n");
    }

    printf("-----------------\n\n");
}

/**
 * print_array is used for debugging purposes.
 *
 * @param array a pointer to a contig. block of integers.
 * @param sizeArray size of the array.
 */
void print_array(int *array, int sizeArray)
{
    for (int i = 0; i < sizeArray; i++)
    {
        printf("%d ", array[i]);
    }

    printf("\n");
}

void signal_sorting_turn()
{
    for (int i = 0; i < NUMBER_OF_THREADS; i++)
    {
        sem_post(&sorting_sem);
    }
}

void wait_for_nonsorting_turn()
{
    for (int i = 0; i < NUMBER_OF_THREADS; i++)
    {
        sem_wait(&nonsorting_sem);
    }
}

/**
 * check_malloc verifies that a block of memory was successfully allocated.
 * Otherwise, it exits the program.
 * @param pointer a pointer to the allocated memory.
 */
void check_malloc(void *pointer)
{
    if (pointer == NULL)
    {
        exit(EXIT_FAILURE);
    }
}

/**
 * free_matrix frees matrix to avoid memory leaks.
 */
void free_matrix()
{
    free(matrix);
}

/**
 * copy_array deep copies content of from_array to to_array.
 * @param from_array array of integers where the contents are copied from.
 * @param to_array array of integers where the contents are copied.
 */
void copy_array(const int *fromArray, int *toArray)
{
    for (int i = 0; i < size; i++)
    {
        toArray[i] = fromArray[i];
    }
}

/**
 * allocate_matrix creates a 2-d array from the
 * contiguous block of integers of the specified length and width.
 * @param array an array of integers
 * @param rowSize the number of rows in 2-d array
 * @param colSize the number of columns in 2-d array
 * @return 2-d array/ pointer to blocks of integers.
 */
int **allocate_matrix(int *array, int rowSize, int colSize)
{
    int **temp = (int **) malloc(rowSize * sizeof(int *));
    check_malloc(temp);

    for (int i = 0; i < rowSize; i++)
    {
        temp[i] = &(array[i * colSize]);
    }

    return temp;
}

/**
 * This function initializes a reusable dissemination barrier.
 * @param dissemination_barrier pointer to the dissemination barrier
 */
void dissemination_barrier_init(dissemination_barrier_t *dissemination_barrier)
{
    dissemination_barrier->arrive = malloc(sizeof(int) * NUMBER_OF_THREADS);
    dissemination_barrier->done = 0;

    for (int i = 0; i < NUMBER_OF_THREADS; i++)
    {
        dissemination_barrier->arrive[i] = 0;
    }
}

/**
 * This is the main logic for the dissemination barrier for thread with
 * workerID.
 * @param dissemination_barrier pointer to the dissemination barrier
 * @param workerID id of the thread/worker
 */
void dissemination_barrier_wait(dissemination_barrier_t *dissemination_barrier, int workerID)
{
    for (int s = 1; s <= number_of_stages; s++)
    {
        // wait for the thread's id to be reset
        while (dissemination_barrier->arrive[workerID] != 0);
        // set the thread's flag
        dissemination_barrier->arrive[workerID] = s;
        // calculate the id of the thread to wait for based on the stage
        int waitFor = (workerID + (1 << (s - 1))) % NUMBER_OF_THREADS;
        while (dissemination_barrier->arrive[waitFor] != s);
        // turn off the flag of the worker it was waiting for
        dissemination_barrier->arrive[waitFor] = 0;
    }
}

void dissemination_barrier_destroy(dissemination_barrier_t *dissemination_barrier)
{
    free((void*)dissemination_barrier->arrive);
}

/**
 * sortColumnSection is the thread's companion function, which
 * represents the section of sorting/work each thread must accomplish.
 * @param arg
 * @return
 */
void *sortColumnSection(void *arg)
{
    while (1)
    {
        sem_wait(&sorting_sem);

        if (task_complete == 1) return NULL;

        worker_data_t *wdt = (worker_data_t *) arg;
        //printf("in %d\n", wdt->id);

        for (int j = wdt->start; j < wdt->end; j++)
        {
            for (int i = 1; i < LENGTH; i++)
            {
                int key = matrix[i][j], k = i - 1;

                while (k >= 0 && matrix[k][j] > key)
                {
                    matrix[k + 1][j] = matrix[k][j];
                    k--;
                }

                matrix[k + 1][j] = key;
            }
        }

        // point in program at which all threads much arrive before any can proceed
        dissemination_barrier_wait(&barrier, wdt->id);
        //printf("out %d\n", wdt->id);
        sem_post(&nonsorting_sem);
    }
}

/**
 * worker_data_init initializes the data that is needed
 * by the thread to know its portion of the work.
 * @param worker_data
 * @param worker_id
 * @param width
 */
void worker_data_init(worker_data_t *worker_data, int worker_id, int width)
{
    worker_data->start = worker_id * work_size;
    worker_data->id = worker_id;

    if (NUMBER_OF_THREADS == (worker_id + 1))
    {
        worker_data->end = width;
    } else
    {
        worker_data->end = worker_data->start + work_size;
    }
}

/**
 * transpose transposes the rxs matrix into an sxr matrix.
 * The it reshapes the matrix back into an rxs matrix.
 * @param length number of rows
 * @param width number of elements in each row
 */
void transpose(int length, int width)
{
    int *auxiliary_array = malloc(size * sizeof(int));
    check_malloc(auxiliary_array);
    int count = 0;

    for (int col = 0; col < width; col++)
    {
        for (int row = 0; row < length; row++)
        {
            auxiliary_array[count] = matrix[row][col];
            count++;
        }
    }

    free_matrix();
    matrix = allocate_matrix(auxiliary_array, length, width);
}

/**
 * reverse_transpose reshapes each set of r/s rows into a single r/s rows
 * into a single r-element row and then transpose the matrix.
 * @param length number of rows
 * @param width number of elements in each row
 */
void reverse_transpose(int length, int width)
{
    int *auxiliary_array = malloc(size * sizeof(int));
    check_malloc(auxiliary_array);
    int **temp = allocate_matrix(auxiliary_array, length, width);
    int a_row = 0, a_col = 0;

    for (int row = 0; row < length; row++)
    {
        for (int col = 0; col < width; col++)
        {
            temp[a_row][a_col] = matrix[row][col];
            a_row++;
            if (a_row == length)
            {
                a_row = 0;
                a_col++;
            }
        }
    }

    free_matrix();
    matrix = allocate_matrix(auxiliary_array, length, width);
}

/**
 * shift_up shifts each column down by r/2 positions, wrapping around into
 * the next column as necessary. The first r/2 entries are filled with MIN.
 * The bottommost r/2 entries are filled with MAX.
 * @param length
 * @param width
 */
void shift_down(int length, int width)
{
    int a_col_size = width + 1, a_row = 0, a_col = 0, shift = length / 2;
    int *auxiliary_array = malloc(a_col_size * length * sizeof(int));
    check_malloc(auxiliary_array);
    int **temp = allocate_matrix(auxiliary_array, length, a_col_size);

    for (int i = 0; i < shift; i++)
    {
        temp[a_row][a_col] = MIN;
        a_row++;
    }

    for (int col = 0; col < width; col++)
    {
        for (int row = 0; row < length; row++)
        {
            temp[a_row][a_col] = matrix[row][col];
            a_row++;
            if (a_row == length)
            {
                a_row = 0;
                a_col++;
            }
        }
    }

    for (int i = 0; i < shift; i++)
    {
        temp[a_row][a_col] = MAX;
        a_row++;
    }

    free_matrix();
    matrix = allocate_matrix(auxiliary_array, length, a_col_size);
}

/**
 * shift_up shifts each column up by r/2 positions,
 * wrapping around into the previous column if needed.
 * @param length
 * @param width
 */
void shift_up(int length, int width)
{
    int col_size = width - 1, shift = length / 2, a_row = 0, a_col = 0;
    int *auxiliary_array = malloc(length * col_size * sizeof(int));
    check_malloc(auxiliary_array);
    int **temp = allocate_matrix(auxiliary_array, length, col_size);

    for (int col = 0; col < width; col++)
    {
        for (int row = 0; row < length; row++)
        {
            if ((col == 0 && row < shift) ||
                (col == col_size && row >= shift))
            {
                continue;
            }

            temp[a_row][a_col] = matrix[row][col];
            a_row++;

            if (a_row == length)
            {
                a_row = 0;
                a_col++;
            }
        }
    }

    free_matrix();
    matrix = allocate_matrix(auxiliary_array, length, col_size);
}

/**
 * copy_matrix_to_array copies the sorted matrix to array A.
 * @param A a pointer to the a block of integers
 * @param length number of rows
 * @param width number of columns
 */
void copy_matrix_to_array(int *A, int lenght, int width)
{
    int count = 0;

    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < lenght; j++)
        {
            A[count] = matrix[j][i];
            count++;
        }
    }
}

/**
 * columnSort implements the parallelized version ofa column sort function that sorts non-negative integers.
 * It follows eight-step algorithm described in the spec.
 * The interface to the function is provided through columnSort.h
 * @param A input array to sort
 * @param numThreads not used for sequential sort
 * @param length number of rows
 * @param width number of columns
 * @param elapsedTime used to pass back the elapsed time to do column-sort.
 */
void columnSort(int *A, int numThreads, int length, int width, double *elapsedTime)
{
    // initialize global variable
    LENGTH = length;
    size = length * width;
    int *B = malloc(size * sizeof(int));
    check_malloc(B);
    copy_array(A, B);
    matrix = allocate_matrix(B, length, width);

    // initialize variables for the parallelization
    NUMBER_OF_THREADS = numThreads;
    if (numThreads > width) NUMBER_OF_THREADS = width;
    work_size = width / numThreads;
    // the number of rounds for the dissemination barrier
    number_of_stages = ceil(log2(numThreads));

    // use start, stop to time the program
    struct timeval start, stop;
    // start timer
    gettimeofday(&start, NULL);

    // initialize semaphores for giving turns to sorting and non-sorting parts.
    sem_init(&nonsorting_sem, 0, 0);
    sem_init(&sorting_sem, 0, NUMBER_OF_THREADS);
    // initialize the dissemination barrier
    dissemination_barrier_init(&barrier);

    // create threads
    pthread_t threadIDs[NUMBER_OF_THREADS];
    worker_data_t workers[NUMBER_OF_THREADS];

    for (int workerId = 0; workerId < NUMBER_OF_THREADS; workerId++)
    {
        worker_data_init(&workers[workerId], workerId, width);
        pthread_create(&threadIDs[workerId], NULL, sortColumnSection, &workers[workerId]);
    }

    // wait for non-sorting turn
    wait_for_nonsorting_turn();
    transpose(length, width);
    // signal that sorting part can enter
    signal_sorting_turn();

    wait_for_nonsorting_turn();
    reverse_transpose(length, width);
    signal_sorting_turn();

    wait_for_nonsorting_turn();
    shift_down(length, width);
    // give the added column to the last worker.
    workers[NUMBER_OF_THREADS - 1].end += 1;
    signal_sorting_turn();

    wait_for_nonsorting_turn();
    shift_up(length, width + 1);
    // used to signal workers to stop
    task_complete = 1;
    signal_sorting_turn();

    // join all threads to the parent thread before the process completes.
    for (int workerId = 0; workerId < NUMBER_OF_THREADS; workerId++)
    {
        pthread_join(threadIDs[workerId], NULL);
    }

    // stop timer
    gettimeofday(&stop, NULL);

    // copy the result of the calculation to the params.
    *elapsedTime = ((stop.tv_sec - start.tv_sec) * 1000000 + (stop.tv_usec - start.tv_usec)) / 1000000.0;
    copy_matrix_to_array(A, length, width);

    // free memory
    sem_destroy(&sorting_sem);
    sem_destroy(&nonsorting_sem);
    dissemination_barrier_destroy(&barrier);
    free_matrix();
}

