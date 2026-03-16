// Accept n inputs from user
// Create a nxn matrix with random ints. Generate random ints within the program
// Take note its time_before
// Transfrom X via a call to the func zsn_optimized(X, n, n)
// Take note of time_after calc
// Obtain elapsed time:: time_elapsed := time_after - time_before
// Show/output elapsed time
#define _POSIX_C_SOURCE 200112L 
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <string.h>
#include <unistd.h>
#include <sched.h>
// #include <linux/time.h>

// function prototypes
double** allocate_matrix(int n);
void free_matrix(double **X, int n);
void generate_rand_matrix(double **X, int n);
void *zsn_optimized(void *arg);
void print_matrix_preview(double **X, int r, int c, int limit);
void *zsn_col_partition(void *arg);
void *zsn_row_partition(void *arg);
double** copy_matrix(double **src, int n);

// thread structure for data computation
// column-partition thread struct
typedef struct {
    double **X;
    int rows;
    int init_col;
    int end_col;
    int thread_id;
} ColThreadData;

// row-partition thread struct
typedef struct {
    double **X;
    int cols;
    int init_row;
    int end_row;
    int thread_id;
} RowThreadData;

int main() {
    int n, num_threads;
    clock_t time_before, time_after;
    double time_elapsed;
    struct timespec start, end;

    printf("Z-Score Normalization for n x n Matrix\n");
    printf("Enter n-size for matrix: ");
    scanf("%d", &n);

    printf("Enter number of threads: ");
    scanf("%d", &num_threads);

    // check if n is positive
    if (n <= 0 || num_threads <= 0) return 1;

    // <--- PERFORM TASK 3 TIMES --->
    // for (int k = 0; k < 3; k++) {
        // random number generator
        srand(time(NULL));
        // call function to allocate matrix using n
        double **X = allocate_matrix(n);
        // generate random nums for matrix nxn
        generate_rand_matrix(X, n);

        // <-- COLUMN PARTITIONING -->    
        // create threads
        double **X_col = copy_matrix(X, n);

        pthread_t col_threads[num_threads];
        ColThreadData col_data[num_threads];

        int cols_per_thread = n / num_threads;

        // take note of time_before
        clock_gettime(CLOCK_MONOTONIC, &start);

        // assign values for each threads and perform zsn_optimized() for each
        for (int i = 0; i < num_threads; i++) {
            col_data[i].X = X_col;
            col_data[i].rows = n;
            col_data[i].init_col = i * cols_per_thread;
            col_data[i].end_col = (i == num_threads - 1) ? n : (i + 1) * cols_per_thread;
            col_data[i].thread_id = i;

            // pthread_create(&threads[i], NULL, zsn_optimized, &data[i]);
            pthread_create(&col_threads[i], NULL, zsn_col_partition, &col_data[i]);
        } 

        for (int i = 0; i < num_threads; i++) {
            pthread_join(col_threads[i], NULL);
        }

        // take note of time_after after running zsn
        clock_gettime(CLOCK_MONOTONIC, &end);

        // print_matrix_preview(X, n, n, 20);

        time_elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;
        printf("Time elapsed for column-wise: %.6f\n", time_elapsed);

        // <--- ROW PARTITIONING --->
        // double **X_row = copy_matrix(X, n);

        // pthread_t row_threads[num_threads];
        // RowThreadData row_data[num_threads];

        // int rows_per_thread = n / num_threads;

        // clock_gettime(CLOCK_MONOTONIC, &start);

        // for (int i = 0; i < num_threads; i++) {
        //     row_data[i].X = X_row;
        //     row_data[i].cols = n;
        //     row_data[i].init_row = i * rows_per_thread;
        //     row_data[i].end_row = (i == num_threads - 1) ? n : (i + 1) * rows_per_thread;
        //     row_data[i].thread_id = i;

        //     pthread_create(&row_threads[i], NULL, zsn_row_partition, &row_data[i]);
        // }

        // for (int i = 0; i < num_threads; i++) {
        //     pthread_join(row_threads[i], NULL);
        // }

        // clock_gettime(CLOCK_MONOTONIC, &end);
        // time_elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;
        // printf("Time elapsed for row-wise: %.6f\n", time_elapsed);

        free_matrix(X, n);
        free_matrix(X_col, n);
        // free_matrix(X_row, n);
    // }

    return 0;    
}

double compute_col_mean(double **X, int r, int c) {
    double sum = 0.0;
    for (int i = 0; i < r; i++) {
        sum += X[i][c];
    }
    return sum / r;
}

double compute_col_std(double **X, int r, int c, double mean) {
    double sum = 0.0;
    for (int i = 0; i < r; i++) {
        double diff = X[i][c] - mean;
        sum += diff * diff;
    }
    return sqrt(sum / r);
}

double** allocate_matrix(int n) {
    double **X = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        X[i] = (double *)malloc(n * sizeof(double));
    }
    return X;
} 

void free_matrix(double **X, int n) {
    for (int i = 0; i < n; i++) {
        free(X[i]);
    }
    free(X);
}

void generate_rand_matrix(double **X, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            X[i][j] = (double)(rand() % 100 + 1);
        }
    }
}

// fully optimized zsn()
// void *zsn_optimized(void *arg) {
//     ThreadData *data = (ThreadData *)arg;

//     for (int i = data->init_col; i < data->end_col; i++) {
//         double sum = 0.0;
//         double sum_sq = 0.0;

//         // compute for statistics in single-pass
//         for (int j = 0; j < data->rows; j++) {
//             double val = data->X[j][i];
//             sum += val;
//             sum += val * val;
//         }

//         double mean = sum / data->rows;
//         double variance = (sum_sq / data->rows) - (mean  * mean);

//         if (variance < 0.0) {
//             variance = 0.0;
//         }

//         double std = sqrt(variance);

//         // Transform with optimized division (multiplication)
//         if (std != 0.0) {
//             double inv_std = 1.0 / std;
//             for (int j = 0; j < data->rows; j++) {
//                 data->X[j][i] = (data->X[j][i] - mean) * inv_std;
//             }
//         } else {
//             for (int j = 0; j < data->rows; j++) {
//                 data->X[j][i] = 0.0;
//             }
//         }
//     }
    
//     return NULL;
// }

void print_matrix_preview(double **X, int r, int c, int limit) {
    int rr = (r < limit) ? r : limit;
    int cc = (c < limit) ? c : limit;

    for (int i = 0; i < rr; i++) {
        for (int j = 0; j < cc; j++)
            printf("%8.4f ", X[i][j]);
        printf("\n");
    }
}

// core-affinity helper
static void pin_thread_core(int thread_id) {
    int num_cores = (int)sysconf(_SC_NPROCESSORS_ONLN);
    int core_id = thread_id % num_cores;

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);

    pthread_t self = pthread_self();
    int rc = pthread_setaffinity_np(self, sizeof(cpu_set_t), &cpuset);
}

// column partitioning
void *zsn_col_partition(void *arg) {
    ColThreadData *data = (ColThreadData *)arg;

    // assign thread to a core
    pin_thread_core(data->thread_id);

    for (int c = data->init_col; c < data->end_col; c++) {
        double sum    = 0.0;
        double sum_sq = 0.0;

        /* Single-pass computation of sum and sum of squares */
        for (int r = 0; r < data->rows; r++) {
            double val = data->X[r][c];
            sum    += val;
            sum_sq += val * val; 
        }

        double mean     = sum / data->rows;
        double variance = (sum_sq / data->rows) - (mean * mean);

        if (variance < 0.0) variance = 0.0;

        double std = sqrt(variance);

        if (std != 0.0) {
            double inv_std = 1.0 / std;
            for (int r = 0; r < data->rows; r++)
                data->X[r][c] = (data->X[r][c] - mean) * inv_std;
        } else {
            for (int r = 0; r < data->rows; r++)
                data->X[r][c] = 0.0;
        }
    }

    return NULL;
}

// row partitioning
void *zsn_row_partition(void *arg) {
    RowThreadData *data = (RowThreadData *)arg;

    // assign thread to a core
    pin_thread_core(data->thread_id);

    for (int r = data->init_row; r < data->end_row; r++) {
        double sum    = 0.0;
        double sum_sq = 0.0;

        /* Single-pass computation of sum and sum of squares */
        for (int c = 0; c < data->cols; c++) {
            double val = data->X[r][c];
            sum    += val;
            sum_sq += val * val;
        }

        double mean     = sum / data->cols;
        double variance = (sum_sq / data->cols) - (mean * mean);

        if (variance < 0.0) variance = 0.0;   /* floating-point guard */

        double std = sqrt(variance);

        if (std != 0.0) {
            double inv_std = 1.0 / std;        /* multiply instead of divide */
            for (int c = 0; c < data->cols; c++)
                data->X[r][c] = (data->X[r][c] - mean) * inv_std;
        } else {
            for (int c = 0; c < data->cols; c++)
                data->X[r][c] = 0.0;
        }
    }

    return NULL;
}

double** copy_matrix(double **src, int n) {
    double **dst = allocate_matrix(n);
    for (int i = 0; i < n; i++)
        memcpy(dst[i], src[i], n * sizeof(double));
    return dst;
}