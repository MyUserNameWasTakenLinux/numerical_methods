# Numerical Methods

A C++ library implementing fundamental numerical methods for solving various mathematical problems.

## Features

1. **Least Squares Solution to an Overdetermined System**  
   Solves overdetermined systems using normal equations to find the best approximation.

2. **QR Factorization Using Householder Transforms**  
   Computes the QR factorization of a given matrix, which is then usedd to solve an overdetermined system.

3. **Eigenvalue and Eigenvector Computation**  
   - **Normalized Power Iteration**: Computes the largest eigenvalue and its corresponding eigenvector.  
   - **Inverse Iteration**: Calculates eigenvalues and eigenvectors with improved precision.

4. **Nonlinear Equation Solvers**  
   - **Bisection Method**: Computes an interval that contains the roots of a continuous function.
   - **Secant Method**: Approximates roots efficiently without requiring derivative information.
  

## Example

**Nonlinear equation solver**
```C++
#include <iostream>
#include <cmath>
#include <utility>
#include "matrix.h"
#include "linear_least_squares.h"
#include "eigen.h"
#include "nonlinear.h"

#define PI 3.14159265

float my_func(float x) {
    return pow(x, 2) - (4 * sin(x));
}

int main(void) {
    auto interval = NONLIN::interval_bisection(1, 3, my_func, 0.1);
    std::cout << "Interval: " << interval.first << " " << interval.second << "\n";
    auto sol = NONLIN::secant_method(1, 3, my_func, 0.0001);
    std::cout << "Secant method solution: " << sol << "\n";
}
```

**Linear Least Squares**
```C++
#include <iostream>
#include <cmath>
#include <utility>
#include "matrix.h"
#include "linear_least_squares.h"
#include "eigen.h"
#include "nonlinear.h"

int main(void) {
    auto A = Matrix(6, 3, {
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
    -1, 1, 0,
    -1, 0, 1,
    0, -1, 1
    });
    
    auto b = Matrix(6, 1, {
    1237,
    1941,
    2417,
    711,
    1177,
    475
    });
    
    auto result = LLS::householder_qr_factorization(A);
    
    for(auto householder_vec : result.first) {
        b = LLS::householder_transform(householder_vec, b);
    }
    // Least Linear Square solution to Ax = b
    auto B = result.second.slice(3, 3, 0, 0);
    auto solution = back_substitution(B, b.slice(3, 1, 0, 0));
    std::cout << "Solution\n" << solution << "\n";
}
```
## Acknowledgments
Algorithms referenced from "Heath, M. T. (2018). Scientific computing : an introductory survey (Second edition, SIAM edition.). Society for Industrial and Applied Mathematics (SIAM, 3600 Market Street, Floor 6, Philadelphia, PA 19104)"
