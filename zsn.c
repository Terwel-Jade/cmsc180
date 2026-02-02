#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function to compute the mean of a column
double compute_column_mean(double **matrix, int rows, int col) {
    double sum = 0.0;
    for (int i = 0; i < rows; i++) {
        sum += matrix[i][col];
    }
    return sum / rows;
}

// Function to compute the standard deviation of a column
double compute_column_std(double **matrix, int rows, int col, double mean) {
    double sum_sq_diff = 0.0;
    for (int i = 0; i < rows; i++) {
        double diff = matrix[i][col] - mean;
        sum_sq_diff += diff * diff;
    }
    return sqrt(sum_sq_diff / rows);
}

// Function to allocate memory for a matrix
double** allocate_matrix(int n) {
    double **matrix = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        matrix[i] = (double *)malloc(n * sizeof(double));
    }
    return matrix;
}

// Function to free matrix memory
void free_matrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// Function to compute Z-Score Normalization
void compute_zsn(double **X, double **T, int n) {
    // For each column j
    for (int j = 0; j < n; j++) {
        // Compute mean (a_j) for column j
        double a_j = compute_column_mean(X, n, j);
        
        // Compute standard deviation (d_j) for column j
        double d_j = compute_column_std(X, n, j, a_j);
        
        // Compute normalized values for column j
        for (int i = 0; i < n; i++) {
            if (d_j != 0.0) {
                T[i][j] = (X[i][j] - a_j) / d_j;
            } else {
                // Handle case where standard deviation is 0 (all values in column are same)
                T[i][j] = 0.0;
            }
        }
    }
}

// Function to print a matrix
void print_matrix(double **matrix, int n, const char *name) {
    printf("\n%s:\n", name);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%8.4f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main() {
    int n;
    
    printf("Z-Score Normalization (ZSN) for n x n Matrix\n");
    printf("=============================================\n\n");
    
    // Input matrix size
    printf("Enter the size of the square matrix (n): ");
    scanf("%d", &n);
    
    if (n <= 0) {
        printf("Error: Matrix size must be positive.\n");
        return 1;
    }
    
    // Allocate memory for matrices X and T
    double **X = allocate_matrix(n);
    double **T = allocate_matrix(n);
    
    // Input matrix X
    printf("\nEnter the elements of the %dx%d matrix X:\n", n, n);
    for (int i = 0; i < n; i++) {
        printf("Row %d: ", i + 1);
        for (int j = 0; j < n; j++) {
            scanf("%lf", &X[i][j]);
        }
    }
    
    // Compute Z-Score Normalization
    compute_zsn(X, T, n);
    
    // Display results
    print_matrix(X, n, "Original Matrix X");
    print_matrix(T, n, "Z-Score Normalized Matrix T");
    
    // Display column statistics
    printf("\nColumn Statistics:\n");
    printf("==================\n");
    for (int j = 0; j < n; j++) {
        double mean = compute_column_mean(X, n, j);
        double std = compute_column_std(X, n, j, mean);
        printf("Column %d: mean = %8.4f, std = %8.4f\n", j + 1, mean, std);
    }
    
    // Free allocated memory
    free_matrix(X, n);
    free_matrix(T, n);
    
    return 0;
}
























#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Function to compute the mean of a column
double compute_column_mean(double **matrix, int rows, int col) {
    double sum = 0.0;
    for (int i = 0; i < rows; i++) {
        sum += matrix[i][col];
    }
    return sum / rows;
}

// Function to compute the standard deviation of a column
double compute_column_std(double **matrix, int rows, int col, double mean) {
    double sum_sq_diff = 0.0;
    for (int i = 0; i < rows; i++) {
        double diff = matrix[i][col] - mean;
        sum_sq_diff += diff * diff;
    }
    return sqrt(sum_sq_diff / rows);
}

// Function to allocate memory for a matrix
double** allocate_matrix(int n) {
    double **matrix = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        matrix[i] = (double *)malloc(n * sizeof(double));
    }
    return matrix;
}

// Function to free matrix memory
void free_matrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// Function to generate random integers for the matrix
void generate_random_matrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Generate random integers between 1 and 1000
            matrix[i][j] = (double)(rand() % 1000 + 1);
        }
    }
}
// Main ZSN function that transforms X in-place (no temporary matrix)
// Following the specification: transform X in place without intermediate computation
// 
// Algorithm:
//   For each column j:
//     1. Compute mean (a_j) and std dev (d_j) for the original column values
//     2. Transform all elements in that column: X[i,j] := (X[i,j] - a_j) / d_j
//   Return transformed matrix X
//
// This approach avoids creating a temporary matrix by processing one column at a time.
// The statistics (mean and std) for each column are computed before any modifications
// to that column, ensuring correct normalization.
void zsn(double **X, int m, int n) {
    // Process each column j
    for (int j = 0; j < n; j++) {
        // Compute mean (a_j) for column j BEFORE any modifications
        double a_j = compute_column_mean(X, m, j);
        
        // Compute standard deviation (d_j) for column j BEFORE any modifications
        double d_j = compute_column_std(X, m, j, a_j);
        
        // Now transform all elements in column j
        for (int i = 0; i < m; i++) {
            // Apply equation 1: X[i,j] := (X[i,j] - a_j) / d_j
            if (d_j != 0.0) {
                X[i][j] = (X[i][j] - a_j) / d_j;
            } else {
                // Handle case where standard deviation is 0
                X[i][j] = 0.0;
            }
        }
    }
    // zsn := X (function returns the transformed matrix X)
}

// Function to print a matrix (only prints if n <= 10 for readability)
void print_matrix(double **matrix, int n, const char *name) {
    if (n <= 10) {
        printf("\n%s:\n", name);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%8.4f ", matrix[i][j]);
            }
            printf("\n");
        }
    } else {
        printf("\n%s: [Matrix too large to display - %dx%d]\n", name, n, n);
        printf("Showing first 5x5 corner:\n");
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                printf("%8.4f ", matrix[i][j]);
            }
            printf("...\n");
        }
        printf("...\n");
    }
}

int main() {
    int n;
    clock_t time_before, time_after;
    double time_elapsed;
    
    printf("Z-Score Normalization (ZSN) for n x n Matrix\n");
    printf("=============================================\n\n");
    
    // Accept n inputs from user
    printf("Enter the size of the square matrix (n): ");
    scanf("%d", &n);
    
    if (n <= 0) {
        printf("Error: Matrix size must be positive.\n");
        return 1;
    }
    
    // Seed the random number generator
    srand(time(NULL));
    
    // Allocate memory for matrix X
    double **X = allocate_matrix(n);
    
    // Create a nxn matrix with random ints
    printf("\nGenerating %dx%d matrix with random integers...\n", n, n);
    generate_random_matrix(X, n);
    
    // Display original matrix (if small enough)
    print_matrix(X, n, "Original Matrix X");
    
    // Take note of time_before
    time_before = clock();
    
    // Transform X via a call to the func zsn(X, n, n)
    zsn(X, n, n);
    
    // Take note of time_after
    time_after = clock();
    
    // Obtain elapsed time
    time_elapsed = ((double)(time_after - time_before)) / CLOCKS_PER_SEC;
    
    // Display normalized matrix (if small enough)
    print_matrix(X, n, "Z-Score Normalized Matrix X");
    
    // Display column statistics (only for small matrices)
    if (n <= 10) {
        printf("\nColumn Statistics (After Normalization):\n");
        printf("========================================\n");
        for (int j = 0; j < n; j++) {
            double mean = compute_column_mean(X, n, j);
            double std = compute_column_std(X, n, j, mean);
            printf("Column %d: mean = %8.4f, std = %8.4f\n", j + 1, mean, std);
        }
    }
    
    // Show/output elapsed time
    printf("\n");
    printf("=============================================\n");
    printf("Matrix size: %d x %d\n", n, n);
    printf("Total elements: %d\n", n * n);
    printf("Time elapsed: %.6f seconds\n", time_elapsed);
    printf("Time elapsed: %.3f milliseconds\n", time_elapsed * 1000);
    printf("=============================================\n");
    
    // Free allocated memory
    free_matrix(X, n);
    
    return 0;
}