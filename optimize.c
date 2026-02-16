#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// ============================================================================
// ORIGINAL TWO-PASS IMPLEMENTATION (Your Code)
// ============================================================================

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

void zsn(double **X, int r, int c) {
    for (int i = 0; i < c; i++) {
        double a = compute_col_mean(X, r, i);
        double b = compute_col_std(X, r, i, a);  // Fixed: was 'c', should be 'i'

        for (int j = 0; j < r; j++) {
            if (b != 0.0) {
                X[j][i] = (X[j][i] - a) / b;
            } else {
                X[j][i] = 0.0;
            }
        }
    }
}

// ============================================================================
// SINGLE-PASS IMPLEMENTATION - METHOD 1: ALGEBRAIC FORMULA
// Reference: Knuth, TAOCP Vol 2, Section 4.2.2
// ============================================================================

/**
 * Compute column mean and standard deviation in a single pass
 * 
 * Mathematical basis:
 *   Variance = E[X²] - E[X]²
 *   Where E[X] = mean, E[X²] = mean of squares
 * 
 * Derivation:
 *   σ² = (1/n)Σ(x_i - μ)²
 *      = (1/n)Σ(x_i² - 2x_iμ + μ²)
 *      = (1/n)Σx_i² - 2μ(1/n)Σx_i + μ²
 *      = (1/n)Σx_i² - 2μ² + μ²
 *      = (1/n)Σx_i² - μ²
 * 
 * References:
 * - Knuth, D.E. (1998). The Art of Computer Programming, Vol 2, Section 4.2.2
 * - Weisstein, E.W. "Sample Variance." MathWorld.
 */
void compute_col_stats_single_pass(double **X, int r, int c, 
                                   double *mean, double *std) {
    double sum = 0.0;      // Σx_i
    double sum_sq = 0.0;   // Σx_i²
    
    // Single pass: accumulate both sums simultaneously
    for (int i = 0; i < r; i++) {
        double val = X[i][c];
        sum += val;
        sum_sq += val * val;
    }
    
    // Calculate mean: μ = (1/n)Σx_i
    *mean = sum / r;
    
    // Calculate variance: σ² = (1/n)Σx_i² - μ²
    double variance = (sum_sq / r) - ((*mean) * (*mean));
    
    // Guard against numerical errors (variance should never be negative)
    // This can happen due to floating-point precision with the algebraic formula
    if (variance < 0.0) {
        variance = 0.0;
    }
    
    // Calculate standard deviation: σ = √(σ²)
    *std = sqrt(variance);
}

void zsn_single_pass(double **X, int r, int c) {
    for (int i = 0; i < c; i++) {
        double a, b;
        
        // Single-pass computation of mean and std
        compute_col_stats_single_pass(X, r, i, &a, &b);

        // Transform column using computed statistics
        for (int j = 0; j < r; j++) {
            if (b != 0.0) {
                X[j][i] = (X[j][i] - a) / b;
            } else {
                X[j][i] = 0.0;
            }
        }
    }
}

// ============================================================================
// SINGLE-PASS IMPLEMENTATION - METHOD 2: WELFORD'S ALGORITHM
// Reference: Welford (1962), Technometrics 4(3):419-420
// ============================================================================

/**
 * Welford's online algorithm for numerically stable variance computation
 * 
 * This method is more numerically stable than the algebraic formula,
 * especially when dealing with large means or small variances.
 * 
 * Algorithm maintains running mean and sum of squared differences:
 *   M_n = M_(n-1) + (x_n - M_(n-1))/n
 *   S_n = S_(n-1) + (x_n - M_(n-1))(x_n - M_n)
 *   Variance = S_n / n
 * 
 * References:
 * - Welford, B.P. (1962). "Note on a method for calculating corrected 
 *   sums of squares and products". Technometrics 4(3): 419-420.
 * - Knuth, D.E. (1998). TAOCP Vol 2, Section 4.2.2, Equation 15
 * - Chan et al. (1983). "Algorithms for Computing the Sample Variance"
 */
void compute_col_stats_welford(double **X, int r, int c,
                               double *mean, double *std) {
    double M = 0.0;  // Running mean
    double S = 0.0;  // Sum of squared differences from current mean
    
    for (int i = 0; i < r; i++) {
        double x = X[i][c];
        
        // Calculate difference from old mean
        double delta = x - M;
        
        // Update mean incrementally
        M += delta / (i + 1);
        
        // Calculate difference from new mean
        double delta2 = x - M;
        
        // Update sum of squared differences
        // This is numerically stable!
        S += delta * delta2;
    }
    
    *mean = M;
    *std = sqrt(S / r);
}

void zsn_welford(double **X, int r, int c) {
    for (int i = 0; i < c; i++) {
        double a, b;
        
        // Welford's algorithm for stable statistics
        compute_col_stats_welford(X, r, i, &a, &b);

        // Transform column
        for (int j = 0; j < r; j++) {
            if (b != 0.0) {
                X[j][i] = (X[j][i] - a) / b;
            } else {
                X[j][i] = 0.0;
            }
        }
    }
}

// ============================================================================
// MAXIMUM OPTIMIZATION: COMBINED + DIVISION TO MULTIPLICATION
// ============================================================================

void zsn_fully_optimized(double **X, int r, int c) {
    for (int col = 0; col < c; col++) {
        double sum = 0.0;
        double sum_sq = 0.0;
        
        // Pass 1: Compute statistics in single pass
        for (int row = 0; row < r; row++) {
            double val = X[row][col];
            sum += val;
            sum_sq += val * val;
        }
        
        double mean = sum / r;
        double variance = (sum_sq / r) - (mean * mean);
        
        if (variance < 0.0) variance = 0.0;
        double std = sqrt(variance);
        
        // Pass 2: Transform with optimized division
        if (std != 0.0) {
            double inv_std = 1.0 / std;  // One division
            for (int row = 0; row < r; row++) {
                X[row][col] = (X[row][col] - mean) * inv_std;  // Multiplication
            }
        } else {
            for (int row = 0; row < r; row++) {
                X[row][col] = 0.0;
            }
        }
    }
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

double** allocate_matrix(int n) {
    double **matrix = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        matrix[i] = (double *)malloc(n * sizeof(double));
    }
    return matrix;
}

void free_matrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void generate_random_matrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = (double)(rand() % 1000 + 1);
        }
    }
}

void copy_matrix(double **dest, double **src, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

// ============================================================================
// MAIN: PERFORMANCE COMPARISON
// ============================================================================

int main() {
    int n;
    clock_t time_before, time_after;
    double time_original, time_single, time_welford, time_optimized;
    
    printf("Single-Pass Statistics Implementation Comparison\n");
    printf("================================================\n\n");
    
    printf("Enter matrix size (n): ");
    if (scanf("%d", &n) != 1 || n <= 0) {
        printf("Error: Invalid matrix size.\n");
        return 1;
    }
    
    srand(time(NULL));
    
    // Allocate four identical matrices
    double **X1 = allocate_matrix(n);
    double **X2 = allocate_matrix(n);
    double **X3 = allocate_matrix(n);
    double **X4 = allocate_matrix(n);
    
    printf("\nGenerating %dx%d random matrix...\n", n, n);
    generate_random_matrix(X1, n);
    copy_matrix(X2, X1, n);
    copy_matrix(X3, X1, n);
    copy_matrix(X4, X1, n);
    
    printf("\n=== PERFORMANCE COMPARISON ===\n");
    
    // Test 1: Original two-pass
    printf("\n1. Original (Two-Pass):\n");
    time_before = clock();
    zsn(X1, n, n);
    time_after = clock();
    time_original = ((double)(time_after - time_before)) / CLOCKS_PER_SEC;
    printf("   Time: %.6f seconds\n", time_original);
    
    // Test 2: Single-pass algebraic
    printf("\n2. Single-Pass (Algebraic Formula):\n");
    time_before = clock();
    zsn_single_pass(X2, n, n);
    time_after = clock();
    time_single = ((double)(time_after - time_before)) / CLOCKS_PER_SEC;
    printf("   Time: %.6f seconds\n", time_single);
    
    // Test 3: Welford's algorithm
    printf("\n3. Single-Pass (Welford's Algorithm):\n");
    time_before = clock();
    zsn_welford(X3, n, n);
    time_after = clock();
    time_welford = ((double)(time_after - time_before)) / CLOCKS_PER_SEC;
    printf("   Time: %.6f seconds\n", time_welford);
    
    // Test 4: Fully optimized
    printf("\n4. Fully Optimized (Single-Pass + Mult):\n");
    time_before = clock();
    zsn_fully_optimized(X4, n, n);
    time_after = clock();
    time_optimized = ((double)(time_after - time_before)) / CLOCKS_PER_SEC;
    printf("   Time: %.6f seconds\n", time_optimized);
    
    // Speedup analysis
    printf("\n=== SPEEDUP ANALYSIS ===\n");
    printf("\nAlgebraic vs Original: %.2fx (%.1f%% faster)\n",
           time_original / time_single,
           (1.0 - time_single/time_original) * 100);
    
    printf("Welford vs Original: %.2fx (%.1f%% faster)\n",
           time_original / time_welford,
           (1.0 - time_welford/time_original) * 100);
    
    printf("Optimized vs Original: %.2fx (%.1f%% faster)\n",
           time_original / time_optimized,
           (1.0 - time_optimized/time_original) * 100);
    
    // Verify correctness
    printf("\n=== CORRECTNESS VERIFICATION ===\n");
    double max_diff = 0.0;
    int all_correct = 1;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double diff1 = fabs(X1[i][j] - X2[i][j]);
            double diff2 = fabs(X1[i][j] - X3[i][j]);
            double diff3 = fabs(X1[i][j] - X4[i][j]);
            
            if (diff1 > max_diff) max_diff = diff1;
            if (diff2 > max_diff) max_diff = diff2;
            if (diff3 > max_diff) max_diff = diff3;
            
            if (diff1 > 1e-10 || diff2 > 1e-10 || diff3 > 1e-10) {
                all_correct = 0;
            }
        }
    }
    
    if (all_correct) {
        printf("✓ All methods produce identical results!\n");
        printf("  Max difference: %.2e (within floating-point precision)\n", max_diff);
    } else {
        printf("✗ Warning: Results differ by more than expected!\n");
        printf("  Max difference: %.2e\n", max_diff);
    }
    
    printf("\n=== ALGORITHM DETAILS ===\n");
    printf("Original:    2 passes per column (mean, then std)\n");
    printf("Algebraic:   1 pass per column (both statistics)\n");
    printf("Welford:     1 pass per column (numerically stable)\n");
    printf("Optimized:   1 pass + multiplication instead of division\n");
    
    printf("\n=== REFERENCES ===\n");
    printf("Algebraic Formula:\n");
    printf("  - Knuth, TAOCP Vol 2 (1998), Section 4.2.2\n");
    printf("  - Weisstein, MathWorld: Sample Variance\n");
    printf("\nWelford's Algorithm:\n");
    printf("  - Welford (1962), Technometrics 4(3):419-420\n");
    printf("  - Chan et al. (1983), The American Statistician 37(3):242-247\n");
    
    // Cleanup
    free_matrix(X1, n);
    free_matrix(X2, n);
    free_matrix(X3, n);
    free_matrix(X4, n);
    
    return 0;
}