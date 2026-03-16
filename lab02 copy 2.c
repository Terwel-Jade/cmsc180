// Accept n inputs from user
// Create a nxn matrix with random ints. Generate random ints within the program
// Take note its time_before
// Transfrom X via a call to the func zsn_optimized(X, n, n)
// Take note of time_after calc
// Obtain elapsed time:: time_elapsed := time_after - time_before
// Show/output elapsed time
#define _POSIX_C_SOURCE 200112L 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
// #include <linux/time.h>

// function prototypes
double** allocate_matrix(int n);
void free_matrix(double **X, int n);
void generate_rand_matrix(double **X, int n);
void *zsn_optimized(void *arg);
void print_matrix_preview(double **X, int r, int c, int limit);

// thread structure for data computation
typedef struct {
    double **X;
    int rows;
    int init_col;
    int end_col;
} ThreadData;

int main() {
    // n is the user input variable
    int n, num_threads;
    // time variables
    clock_t time_before, time_after;
    double time_elapsed;
    // time vars for threaded processes
    struct timespec start, end;

    printf("Z-Score Normalization for n x n Matrix\n");
    printf("Enter n-size for matrix: ");
    scanf("%d", &n);

    printf("Enter number of threads: ");
    scanf("%d", &num_threads);

    // check if n is positive
    if (n < 0 || num_threads <= 0) return 1;

    // random number generator
    srand(time(NULL));

    // call function to allocate matrix using n
    double **X = allocate_matrix(n);

    // generate random nums for matrix nxn
    generate_rand_matrix(X, n);

    // create threads
    pthread_t threads[num_threads];
    ThreadData data[num_threads];

    int cols_per_thread = n / num_threads;

    // print_matrix_preview(X, n, n, 20);
    // printf("\n");

    // take note of time_before
    clock_gettime(CLOCK_MONOTONIC, &start);

    // assign values for each threads and perform zsn_optimized() for each
    for (int i = 0; i < num_threads; i++) {
        data[i].X = X;
        data[i].rows = n;
        data[i].init_col = i * cols_per_thread;
        data[i].end_col = (i == num_threads - 1) ? n : (i + 1) * cols_per_thread;

        // pthread_create(&threads[i], NULL, zsn_optimized, &data[i]);
        pthread_create(&threads[i], NULL, zsn_optimized, &data[i]);
    } 

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    // take note of time_after after running zsn
    clock_gettime(CLOCK_MONOTONIC, &end);

    // print_matrix_preview(X, n, n, 20);

    // time_elapsed = ((double)(time_after - time_before)) / CLOCKS_PER_SEC;
    time_elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;
    printf("Time elapsed: %.6f\n", time_elapsed);

    free_matrix(X, n);

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
void *zsn_optimized(void *arg) {
    ThreadData *data = (ThreadData *)arg;

    for (int i = data->init_col; i < data->end_col; i++) {
        double sum = 0.0;
        double sum_sq = 0.0;

        // compute for statistics in single-pass
        for (int j = 0; j < data->rows; j++) {
            double val = data->X[j][i];
            sum += val;
            sum += val * val;
        }

        double mean = sum / data->rows;
        double variance = (sum_sq / data->rows) - (mean  * mean);

        if (variance < 0.0) {
            variance = 0.0;
        }

        double std = sqrt(variance);

        // Transform with optimized division (multiplication)
        if (std != 0.0) {
            double inv_std = 1.0 / std;
            for (int j = 0; j < data->rows; j++) {
                data->X[j][i] = (data->X[j][i] - mean) * inv_std;
            }
        } else {
            for (int j = 0; j < data->rows; j++) {
                data->X[j][i] = 0.0;
            }
        }
    }
    
    return NULL;
}

void print_matrix_preview(double **X, int r, int c, int limit) {
    int rr = (r < limit) ? r : limit;
    int cc = (c < limit) ? c : limit;

    for (int i = 0; i < rr; i++) {
        for (int j = 0; j < cc; j++)
            printf("%8.4f ", X[i][j]);
        printf("\n");
    }
}