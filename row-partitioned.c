/*
 * Z-Score Normalization (ZSN) for n x n Matrix
 *
 * Features:
 *   - Column-partitioning  : each thread normalizes a subset of columns
 *   - Row-partitioning     : each thread normalizes a subset of rows
 *   - Core affinity        : each thread is pinned to a dedicated CPU core
 *                            (step 7 from the lab spec — very important)
 *
 * Build:
 *   gcc -O2 -o zsn_threaded zsn_threaded.c -lm -lpthread
 *
 * Usage:
 *   ./zsn_threaded          (interactive prompts)
 *   ./zsn_threaded <n> <t>  (command-line / data-stream mode)
 */

#define _GNU_SOURCE          /* needed for pthread_setaffinity_np & CPU_* */
#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <sched.h>          /* cpu_set_t, CPU_SET, CPU_ZERO              */

/* ------------------------------------------------------------------ */
/*  Thread argument structs                                             */
/* ------------------------------------------------------------------ */

/* Column-partition: thread owns columns [init_col, end_col)           */
typedef struct {
    double **X;
    int      rows;      /* == n                                        */
    int      init_col;
    int      end_col;
    int      thread_id; /* used for core affinity                      */
} ColThreadData;

/* Row-partition: thread owns rows [init_row, end_row)                 */
typedef struct {
    double **X;
    int      cols;      /* == n                                        */
    int      init_row;
    int      end_row;
    int      thread_id; /* used for core affinity                      */
} RowThreadData;

/* ------------------------------------------------------------------ */
/*  Prototypes                                                         */
/* ------------------------------------------------------------------ */
double** allocate_matrix(int n);
void     free_matrix(double **X, int n);
void     generate_rand_matrix(double **X, int n);

/* Core-affinity helper */
static void pin_thread_to_core(int thread_id);

/* Column-partition ZSN (original, now with affinity) */
void    *zsn_col_partition(void *arg);

/* Row-partition ZSN (new) */
void    *zsn_row_partition(void *arg);

/* Utility */
void     print_matrix_preview(double **X, int r, int c, int limit);
double** copy_matrix(double **src, int n);

/* ------------------------------------------------------------------ */
/*  main                                                               */
/* ------------------------------------------------------------------ */
int main(int argc, char *argv[]) {
    int n, num_threads;
    struct timespec start, end;
    double time_elapsed;

    /* ----- Read n and t (command-line or interactive) -------------- */
    if (argc == 3) {
        n           = atoi(argv[1]);
        num_threads = atoi(argv[2]);
    } else {
        printf("Z-Score Normalization for n x n Matrix\n");
        printf("Enter n (matrix size)   : ");
        scanf("%d", &n);
        printf("Enter t (num threads)   : ");
        scanf("%d", &num_threads);
    }

    if (n <= 0 || num_threads <= 0) {
        fprintf(stderr, "Error: n and t must be positive integers.\n");
        return 1;
    }

    /* Warn if n is not >> t (lab spec reminder) */
    if (n <= num_threads) {
        fprintf(stderr,
            "Warning: n (%d) should be >> t (%d) for meaningful partitioning.\n",
            n, num_threads);
    }

    /* ----- Seed RNG and build matrix ------------------------------- */
    srand((unsigned)time(NULL));
    double **X_original = allocate_matrix(n);
    generate_rand_matrix(X_original, n);

    printf("\nMatrix size : %d x %d\n", n, n);
    printf("Threads     : %d\n\n", num_threads);

    /* ================================================================
     * 1.  COLUMN-PARTITION  (original algorithm, now with affinity)
     * ================================================================ */
    double **X_col = copy_matrix(X_original, n);

    pthread_t    col_threads[num_threads];
    ColThreadData col_data[num_threads];

    int cols_per_thread = n / num_threads;

    clock_gettime(CLOCK_MONOTONIC, &start);

    for (int i = 0; i < num_threads; i++) {
        col_data[i].X         = X_col;
        col_data[i].rows      = n;
        col_data[i].init_col  = i * cols_per_thread;
        col_data[i].end_col   = (i == num_threads - 1)
                                  ? n
                                  : (i + 1) * cols_per_thread;
        col_data[i].thread_id = i;

        pthread_create(&col_threads[i], NULL, zsn_col_partition, &col_data[i]);
    }

    for (int i = 0; i < num_threads; i++)
        pthread_join(col_threads[i], NULL);

    clock_gettime(CLOCK_MONOTONIC, &end);
    time_elapsed = (end.tv_sec  - start.tv_sec)
                 + (end.tv_nsec - start.tv_nsec) * 1e-9;

    printf("=== Column-partition ZSN ===\n");
    printf("Time elapsed : %.9f s\n\n", time_elapsed);

    /* ================================================================
     * 2.  ROW-PARTITION  (new)
     * ================================================================ */
    double **X_row = copy_matrix(X_original, n);

    pthread_t    row_threads[num_threads];
    RowThreadData row_data[num_threads];

    int rows_per_thread = n / num_threads;

    clock_gettime(CLOCK_MONOTONIC, &start);

    for (int i = 0; i < num_threads; i++) {
        row_data[i].X         = X_row;
        row_data[i].cols      = n;
        row_data[i].init_row  = i * rows_per_thread;
        row_data[i].end_row   = (i == num_threads - 1)
                                  ? n
                                  : (i + 1) * rows_per_thread;
        row_data[i].thread_id = i;

        pthread_create(&row_threads[i], NULL, zsn_row_partition, &row_data[i]);
    }

    for (int i = 0; i < num_threads; i++)
        pthread_join(row_threads[i], NULL);

    clock_gettime(CLOCK_MONOTONIC, &end);
    time_elapsed = (end.tv_sec  - start.tv_sec)
                 + (end.tv_nsec - start.tv_nsec) * 1e-9;

    printf("=== Row-partition ZSN ===\n");
    printf("Time elapsed : %.9f s\n\n", time_elapsed);

    /* ----- Optional: preview first few cells ----------------------- */
    if (n <= 100) {
        printf("--- Column-partition result (preview, up to 8x8) ---\n");
        print_matrix_preview(X_col, n, n, 8);
        printf("\n--- Row-partition result (preview, up to 8x8) ---\n");
        print_matrix_preview(X_row, n, n, 8);
        printf("\n");
    }

    /* ----- Cleanup ------------------------------------------------- */
    free_matrix(X_original, n);
    free_matrix(X_col,      n);
    free_matrix(X_row,      n);

    return 0;
}

/* ------------------------------------------------------------------ */
/*  Core-affinity helper                                               */
/*                                                                     */
/*  Pins the calling thread to CPU core (thread_id % available_cores) */
/*  following the lab spec: "Assign the thread to a core"             */
/* ------------------------------------------------------------------ */
static void pin_thread_to_core(int thread_id) {
    int num_cores = (int)sysconf(_SC_NPROCESSORS_ONLN);
    int core_id   = thread_id % num_cores;

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);

    pthread_t self = pthread_self();
    int rc = pthread_setaffinity_np(self, sizeof(cpu_set_t), &cpuset);
    if (rc != 0) {
        /* Non-fatal: print a warning but continue */
        fprintf(stderr,
            "[Thread %d] Warning: could not set affinity to core %d (rc=%d)\n",
            thread_id, core_id, rc);
    }
}

/* ------------------------------------------------------------------ */
/*  Column-partition ZSN  (optimized, with core affinity)             */
/*                                                                     */
/*  Each thread normalizes columns [init_col, end_col).               */
/*  Normalization is column-wise:                                      */
/*    mean  = avg of all rows in that column                           */
/*    std   = population std-dev of that column                        */
/*    X[r][c] = (X[r][c] - mean) / std                                */
/* ------------------------------------------------------------------ */
void *zsn_col_partition(void *arg) {
    ColThreadData *data = (ColThreadData *)arg;

    /* --- STEP 7: Assign thread to a core (very important) --------- */
    pin_thread_to_core(data->thread_id);

    for (int c = data->init_col; c < data->end_col; c++) {
        double sum    = 0.0;
        double sum_sq = 0.0;

        /* Single-pass computation of sum and sum of squares */
        for (int r = 0; r < data->rows; r++) {
            double val = data->X[r][c];
            sum    += val;
            sum_sq += val * val;   /* BUG FIX: was "sum += val*val" in original */
        }

        double mean     = sum / data->rows;
        double variance = (sum_sq / data->rows) - (mean * mean);

        if (variance < 0.0) variance = 0.0;   /* floating-point guard */

        double std = sqrt(variance);

        if (std != 0.0) {
            double inv_std = 1.0 / std;        /* multiply instead of divide */
            for (int r = 0; r < data->rows; r++)
                data->X[r][c] = (data->X[r][c] - mean) * inv_std;
        } else {
            for (int r = 0; r < data->rows; r++)
                data->X[r][c] = 0.0;
        }
    }

    return NULL;
}

/* ------------------------------------------------------------------ */
/*  Row-partition ZSN  (new, with core affinity)                      */
/*                                                                     */
/*  Each thread normalizes rows [init_row, end_row).                  */
/*  Normalization is row-wise:                                         */
/*    mean  = avg of all columns in that row                           */
/*    std   = population std-dev of that row                           */
/*    X[r][c] = (X[r][c] - mean) / std                                */
/* ------------------------------------------------------------------ */
void *zsn_row_partition(void *arg) {
    RowThreadData *data = (RowThreadData *)arg;

    /* --- STEP 7: Assign thread to a core (very important) --------- */
    pin_thread_to_core(data->thread_id);

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

/* ------------------------------------------------------------------ */
/*  Matrix helpers                                                     */
/* ------------------------------------------------------------------ */
double** allocate_matrix(int n) {
    double **X = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
        X[i] = (double *)malloc(n * sizeof(double));
    return X;
}

void free_matrix(double **X, int n) {
    for (int i = 0; i < n; i++)
        free(X[i]);
    free(X);
}

void generate_rand_matrix(double **X, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            X[i][j] = (double)(rand() % 100 + 1);
}

/* Deep copy: allocate a new n×n matrix and copy src into it */
double** copy_matrix(double **src, int n) {
    double **dst = allocate_matrix(n);
    for (int i = 0; i < n; i++)
        memcpy(dst[i], src[i], n * sizeof(double));
    return dst;
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