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
























#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

// Function prototypes
double** allocate_matrix(int n);
void free_matrix(double **X, int n);
void generate_rand_matrix(double **X, int n);
void copy_matrix(double **dest, double **src, int n);
void print_matrix_preview(double **X, int r, int c, int limit);
int verify_results(double **X1, double **X2, int n, double tolerance);

// Single-threaded fully optimized version
void zsn_fully_optimized(double **X, int r, int c);

// Multi-threaded fully optimized version
void *zsn_thread_optimized(void *arg);
void zsn_threaded(double **X, int r, int c, int num_threads);

// Thread structure for data computation
typedef struct {
    double **X;
    int rows;
    int start_col;
    int end_col;
} ThreadData;

int main() {
    int n, num_threads;
    clock_t time_before, time_after;
    double time_single, time_threaded;

    printf("Threaded Fully-Optimized ZSN Implementation\n");
    printf("============================================\n\n");
    
    printf("Enter n-size for matrix: ");
    scanf("%d", &n);

    printf("Enter number of threads: ");
    scanf("%d", &num_threads);

    if (n <= 0 || num_threads <= 0) {
        printf("Error: n and num_threads must be positive\n");
        return 1;
    }

    // Seed random number generator
    srand(time(NULL));

    // Allocate two matrices for comparison
    printf("\nAllocating matrices...\n");
    double **X_single = allocate_matrix(n);
    double **X_threaded = allocate_matrix(n);

    // Generate random matrix
    printf("Generating %dx%d random matrix...\n", n, n);
    generate_rand_matrix(X_single, n);
    copy_matrix(X_threaded, X_single, n);

    printf("\nOriginal matrix preview (first 10x10):\n");
    print_matrix_preview(X_single, n, n, 10);
    printf("\n");

    // ========================================================================
    // TEST 1: SINGLE-THREADED FULLY OPTIMIZED
    // ========================================================================
    printf("=== TEST 1: SINGLE-THREADED (Fully Optimized) ===\n");
    
    time_before = clock();
    zsn_fully_optimized(X_single, n, n);
    time_after = clock();
    time_single = ((double)(time_after - time_before)) / CLOCKS_PER_SEC;
    
    printf("Time: %.6f seconds (%.3f ms)\n", time_single, time_single * 1000);
    printf("\nResult preview:\n");
    print_matrix_preview(X_single, n, n, 10);
    printf("\n");

    // ========================================================================
    // TEST 2: MULTI-THREADED FULLY OPTIMIZED
    // ========================================================================
    printf("=== TEST 2: MULTI-THREADED (Fully Optimized) ===\n");
    printf("Threads: %d\n", num_threads);
    
    time_before = clock();
    zsn_threaded(X_threaded, n, n, num_threads);
    time_after = clock();
    time_threaded = ((double)(time_after - time_before)) / CLOCKS_PER_SEC;
    
    printf("Time: %.6f seconds (%.3f ms)\n", time_threaded, time_threaded * 1000);
    printf("\nResult preview:\n");
    print_matrix_preview(X_threaded, n, n, 10);
    printf("\n");

    // ========================================================================
    // VERIFICATION
    // ========================================================================
    printf("=== VERIFICATION ===\n");
    
    if (verify_results(X_single, X_threaded, n, 1e-10)) {
        printf("✓ SUCCESS: Results are identical!\n");
    } else {
        printf("✗ ERROR: Results differ!\n");
    }

    // ========================================================================
    // PERFORMANCE ANALYSIS
    // ========================================================================
    printf("\n=== PERFORMANCE ANALYSIS ===\n");
    printf("Single-threaded: %.6f seconds\n", time_single);
    printf("Multi-threaded:  %.6f seconds\n", time_threaded);
    
    if (time_threaded < time_single) {
        double speedup = time_single / time_threaded;
        double efficiency = (speedup / num_threads) * 100;
        
        printf("\nSpeedup: %.2fx faster\n", speedup);
        printf("Parallel efficiency: %.1f%%\n", efficiency);
        printf("Time saved: %.6f seconds (%.1f%% improvement)\n",
               time_single - time_threaded,
               (1.0 - time_threaded/time_single) * 100);
    } else {
        printf("\nNo speedup (threading overhead > benefit)\n");
        printf("Try larger n for threading benefits\n");
    }

    printf("\n=== SUMMARY ===\n");
    printf("Matrix size: %dx%d (%d elements)\n", n, n, n*n);
    printf("Columns per thread: ~%d\n", n / num_threads);
    printf("Algorithm: Algebraic formula (E[X²] - E[X]²)\n");
    printf("Optimization: Division to multiplication\n");

    // Cleanup
    free_matrix(X_single, n);
    free_matrix(X_threaded, n);

    return 0;
}

// ============================================================================
// SINGLE-THREADED FULLY OPTIMIZED VERSION
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
// MULTI-THREADED FULLY OPTIMIZED VERSION
// ============================================================================

/**
 * Thread worker function - processes a range of columns
 * Uses algebraic formula for single-pass statistics
 */
void *zsn_thread_optimized(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    
    // Each thread processes its assigned columns
    for (int col = data->start_col; col < data->end_col; col++) {
        double sum = 0.0;
        double sum_sq = 0.0;
        
        // Pass 1: Compute statistics in single pass
        // Using algebraic formula: Var(X) = E[X²] - E[X]²
        for (int row = 0; row < data->rows; row++) {
            double val = data->X[row][col];
            sum += val;
            sum_sq += val * val;
        }
        
        // Calculate mean
        double mean = sum / data->rows;
        
        // Calculate variance using computational formula
        double variance = (sum_sq / data->rows) - (mean * mean);
        
        // Guard against numerical errors
        if (variance < 0.0) {
            variance = 0.0;
        }
        
        // Calculate standard deviation
        double std = sqrt(variance);
        
        // Pass 2: Transform column with optimized division-to-multiplication
        if (std != 0.0) {
            double inv_std = 1.0 / std;  // One division per column
            
            for (int row = 0; row < data->rows; row++) {
                data->X[row][col] = (data->X[row][col] - mean) * inv_std;  // Multiply
            }
        } else {
            // Handle zero variance case
            for (int row = 0; row < data->rows; row++) {
                data->X[row][col] = 0.0;
            }
        }
    }
    
    return NULL;
}

/**
 * Main threaded ZSN function
 * Distributes columns across threads
 */
void zsn_threaded(double **X, int r, int c, int num_threads) {
    pthread_t threads[num_threads];
    ThreadData thread_data[num_threads];
    
    int cols_per_thread = c / num_threads;
    
    // Create threads
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].X = X;
        thread_data[i].rows = r;
        thread_data[i].start_col = i * cols_per_thread;
        
        // Last thread handles remaining columns
        if (i == num_threads - 1) {
            thread_data[i].end_col = c;
        } else {
            thread_data[i].end_col = (i + 1) * cols_per_thread;
        }
        
        pthread_create(&threads[i], NULL, zsn_thread_optimized, &thread_data[i]);
    }
    
    // Wait for all threads to complete
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

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
            X[i][j] = (double)(rand() % 1000 + 1);
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

void print_matrix_preview(double **X, int r, int c, int limit) {
    int rr = (r < limit) ? r : limit;
    int cc = (c < limit) ? c : limit;

    for (int i = 0; i < rr; i++) {
        for (int j = 0; j < cc; j++) {
            printf("%8.4f ", X[i][j]);
        }
        printf("\n");
    }
}

int verify_results(double **X1, double **X2, int n, double tolerance) {
    double max_diff = 0.0;
    int errors = 0;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double diff = fabs(X1[i][j] - X2[i][j]);
            if (diff > max_diff) {
                max_diff = diff;
            }
            if (diff > tolerance) {
                errors++;
                if (errors <= 3) {
                    printf("  Diff at [%d][%d]: %.10f vs %.10f (%.2e)\n",
                           i, j, X1[i][j], X2[i][j], diff);
                }
            }
        }
    }
    
    printf("  Max difference: %.2e\n", max_diff);
    
    if (errors > 0) {
        printf("  Total errors: %d (> tolerance)\n", errors);
        return 0;
    }
    
    return 1;
}