# Mathematical Derivation: Single-Pass Variance Formula

## Visual Proof of E[X²] - E[X]²

### Traditional Formula (Two-Pass Required)
```
σ² = (1/n) Σ(x_i - μ)²

where μ = (1/n) Σx_i

This requires:
  Pass 1: Calculate μ
  Pass 2: Calculate Σ(x_i - μ)²
```

### Computational Formula (Single-Pass)
```
σ² = (1/n) Σx_i² - μ²

This requires:
  Pass 1: Calculate both Σx_i and Σx_i² simultaneously!
```

---

## Step-by-Step Algebraic Proof

### Starting Point
```
σ² = (1/n) Σ(x_i - μ)²
```

### Step 1: Expand the square
```
σ² = (1/n) Σ(x_i² - 2x_iμ + μ²)
```

### Step 2: Distribute the summation
```
σ² = (1/n)[Σx_i² - Σ(2x_iμ) + Σμ²]
```

### Step 3: Factor out constants
```
σ² = (1/n)[Σx_i² - 2μΣx_i + nμ²]
```
Note: 
- μ is constant, so Σ(2x_iμ) = 2μΣx_i
- Σμ² = nμ² (summing μ² exactly n times)

### Step 4: Substitute μ = (1/n)Σx_i
```
σ² = (1/n)[Σx_i² - 2μ·(nμ) + nμ²]

σ² = (1/n)[Σx_i² - 2nμ² + nμ²]
```

### Step 5: Simplify
```
σ² = (1/n)[Σx_i² - nμ²]

σ² = (1/n)Σx_i² - μ²
```

### Final Form
```
σ² = E[X²] - (E[X])²

where:
  E[X²] = (1/n)Σx_i²  (mean of squares)
  E[X]  = (1/n)Σx_i   (mean)
  
Standard deviation:
  σ = √σ²
```

---

## Numerical Example

### Dataset: {2, 4, 6, 8, 10}

#### Method 1: Traditional (Two-Pass)

**Pass 1: Calculate mean**
```
μ = (2 + 4 + 6 + 8 + 10) / 5
μ = 30 / 5 = 6
```

**Pass 2: Calculate variance**
```
σ² = [(2-6)² + (4-6)² + (6-6)² + (8-6)² + (10-6)²] / 5
σ² = [16 + 4 + 0 + 4 + 16] / 5
σ² = 40 / 5 = 8
σ = √8 ≈ 2.828
```

#### Method 2: Computational (Single-Pass)

**Single Pass: Calculate both sums**
```
Σx_i  = 2 + 4 + 6 + 8 + 10 = 30
Σx_i² = 4 + 16 + 36 + 64 + 100 = 220

μ = 30 / 5 = 6
```

**Calculate variance**
```
σ² = (Σx_i² / n) - μ²
σ² = (220 / 5) - 6²
σ² = 44 - 36 = 8
σ = √8 ≈ 2.828
```

**✓ Same result, but only one pass through the data!**

---

## Visual Representation

```
Traditional Method (Two Passes):
=====================================

Data:     [2] [4] [6] [8] [10]
          ↓   ↓   ↓   ↓   ↓
Pass 1:   Sum all → μ = 6
          
Data:     [2] [4] [6] [8] [10]
          ↓   ↓   ↓   ↓   ↓
Pass 2:   (x-μ)² for each → σ² = 8


Computational Method (Single Pass):
=====================================

Data:     [2] [4] [6] [8] [10]
          ↓   ↓   ↓   ↓   ↓
Pass 1:   Σx = 30  (accumulate x)
          Σx² = 220 (accumulate x² simultaneously!)
          
Then:     μ = 30/5 = 6
          σ² = 220/5 - 6² = 44 - 36 = 8
```

---

## Code Comparison

### Traditional (Your Original Code)
```c
// Pass 1: Mean
double sum = 0.0;
for (int i = 0; i < r; i++) {
    sum += X[i][c];
}
double mean = sum / r;

// Pass 2: Variance
double sum_sq_diff = 0.0;
for (int i = 0; i < r; i++) {
    double diff = X[i][c] - mean;
    sum_sq_diff += diff * diff;
}
double variance = sum_sq_diff / r;
double std = sqrt(variance);
```

### Computational (Single-Pass)
```c
// Single pass for both!
double sum = 0.0;
double sum_sq = 0.0;

for (int i = 0; i < r; i++) {
    double val = X[i][c];
    sum += val;      // Σx
    sum_sq += val * val;  // Σx²
}

double mean = sum / r;
double variance = (sum_sq / r) - (mean * mean);
double std = sqrt(variance);
```

---

## Complexity Analysis

### Time Complexity
```
Traditional:  O(2n) per column  (two passes)
Computational: O(n) per column  (one pass)

For nxn matrix:
Traditional:  O(2n²) total
Computational: O(n²) total

Speedup: 2x theoretical, ~1.5x practical (due to cache)
```

### Space Complexity
```
Both methods: O(1) extra space (just scalar variables)
```

---

## Important Notes

### Numerical Stability Warning

The formula `σ² = E[X²] - E[X]²` can suffer from **catastrophic cancellation** when:
- Mean is very large (e.g., 10⁹)
- Variance is very small relative to mean

**Example of problematic case:**
```
Data: [1000000000, 1000000001, 1000000002]

E[X] = 1000000001
E[X²] = 1000000002000000001
E[X²] - (E[X])² might lose precision!
```

**Solution:** Use Welford's algorithm for such cases (see references)

### When to use each method

**Computational Formula (E[X²] - E[X]²):**
- ✓ General purpose (works for 99% of cases)
- ✓ Fastest implementation
- ✓ Simple to understand
- ✗ Can have precision issues with extreme values

**Welford's Algorithm:**
- ✓ Numerically stable for all cases
- ✓ Handles streaming data
- ✓ Better for extreme values
- ✗ Slightly more complex
- ✗ Marginally slower (very small difference)

---

## References for Mathematical Proof

1. **Weisstein, Eric W.** "Sample Variance." From MathWorld--A Wolfram Web Resource.
   https://mathworld.wolfram.com/SampleVariance.html
   - Clear mathematical derivation

2. **Wikipedia** - "Algebraic formula for the variance"
   https://en.wikipedia.org/wiki/Variance#Definition
   - Includes proof of equivalence

3. **Knuth, Donald E.** (1998). The Art of Computer Programming, Vol 2, Section 4.2.2
   - Rigorous treatment of numerical algorithms

4. **Rice, John A.** (2006). Mathematical Statistics and Data Analysis (3rd ed.)
   - Standard textbook proof

5. **MIT OpenCourseWare** - 18.05 Introduction to Probability and Statistics
   - Educational derivation with examples