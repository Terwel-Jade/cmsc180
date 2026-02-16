#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

/* ===================== Function Prototypes ===================== */
double compute_col_mean(double **X, int r, int c);
double compute_col_std(double **X, int r, int c, double mean);

double **allocate_matrix(int n);
void free_matrix(double **X, int n);
void generate_rand_matrix(double **X, int n);

void print_matrix_preview(double **X, int r, int c, int limit);

/* ===================== Thread Data ===================== */
typedef struct {
    double **X;
    int rows;
    int start_col;
    int end_col;
} ThreadData;

/* ===================== Column Functions ===================== */
double compute_col_mean(double **X, int r, int c) {
    double sum = 0.0;
    for (int i = 0; i < r; i++)
        sum += X[i][c];
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

/* ===================== Thread Worker ===================== */
void *zsn_thread(void *arg) {
    ThreadData *data = (ThreadData *)arg;

    for (int i = data->start_col; i < data->end_col; i++) {
        double mean = compute_col_mean(data->X, data->rows, i);
        double std  = compute_col_std(data->X, data->rows, i, mean);

        for (int j = 0; j < data->rows; j++) {
            if (std != 0.0)
                data->X[j][i] = (data->X[j][i] - mean) / std;
            else
                data->X[j][i] = 0.0;
        }
    }
    return NULL;
}

/* ===================== Matrix Utilities ===================== */
double **allocate_matrix(int n) {
    double **X = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
        X[i] = malloc(n * sizeof(double));
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
            X[i][j] = (double)(rand() % 2000 + 1);
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

/* ===================== Main ===================== */
int main() {
    int n, num_threads;
    clock_t start, end;

    printf("Enter matrix size n: ");
    scanf("%d", &n);

    printf("Enter number of threads: ");
    scanf("%d", &num_threads);

    if (n <= 0 || num_threads <= 0) return 1;

    srand(time(NULL));

    double **X = allocate_matrix(n);
    generate_rand_matrix(X, n);

    printf("\nBefore ZSN (preview):\n");
    print_matrix_preview(X, n, n, 5);

    pthread_t threads[num_threads];
    ThreadData data[num_threads];

    int cols_per_thread = n / num_threads;

    start = clock();

    for (int i = 0; i < num_threads; i++) {
        data[i].X = X;
        data[i].rows = n;
        data[i].start_col = i * cols_per_thread;
        data[i].end_col = (i == num_threads - 1)
                            ? n
                            : (i + 1) * cols_per_thread;

        pthread_create(&threads[i], NULL, zsn_thread, &data[i]);
    }

    for (int i = 0; i < num_threads; i++)
        pthread_join(threads[i], NULL);

    end = clock();

    printf("\nAfter ZSN (preview):\n");
    print_matrix_preview(X, n, n, 5);

    printf("\nTime elapsed: %.6f seconds\n",
           (double)(end - start) / CLOCKS_PER_SEC);

    free_matrix(X, n);
    return 0;
}
