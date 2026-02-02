// Accept n inputs from user
// Create a nxn matrix with random ints. Generate random ints within the program
// Take note its time_before
// Transfrom X via a call to the func zsn(X, n, n)
// Take note of time_after calc
// Obtain elapsed time:: time_elapsed := time_after - time_before
// Show/output elapsed time

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// function prototypes
double compute_col_mean(double **X, int r, int c);
double compute_col_std(double **X, int r, int c, double mean);
double** allocate_matrix(int n);
void free_matrix(double **X, int n);
void generate_rand_matrix(double **X, int n);
void zsn(double **X, int r, int c);

int main() {
    // n is the user input variable
    int n;
    // time variables
    clock_t time_before, time_after;
    double time_elapsed;

    printf("Z-Score Normalization for n x n Matrix\n");
    printf("Enter n-size for matrix: ");
    scanf("%d", &n);

    // check if n is positive
    if (n < 0) return 1;

    // random number generator
    srand(time(NULL));

    // call function to allocate matrix using n
    double **X = allocate_matrix(n);

    // generate random nums for matrix nxn
    generate_rand_matrix(X, n);

    // take note of time_before
    time_before = clock();

    // do the calculations; Transform X via zsn();
    zsn(X, n, n);

    // take note of time_after after running zsn
    time_after = clock();

    time_elapsed = ((double)(time_after - time_before)) / CLOCKS_PER_SEC;
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
    for (int i = 0; i < r, i++) {
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
            X[i][j] = (double)(rand() % 2000 + 1);
        }
    }
}

void zsn(double **X, int r, int c) {
    for (int i = 0; i < c; i++) {
        double a = compute_col_mean(X, r, i);
        double b = compute_col_std(X, c, a);

        for (int j = 0; j < r; j++) {
            if (b != 0.0) {
                X[j][i] = (X[j][i] - a) / b;
            } else {
                X[j][i] = 0.0;
            }
        }
    }
}