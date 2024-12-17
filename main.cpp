#include <iostream>
#include <cmath>
#include <utility>
#include "matrix.h"
#include "linear_least_squares.h"
#include "eigen.h"
#include "nonlinear.h"

float f(float x) {
    return (4 / (1 + pow(x, 2)));
} 

int main(void) {
    auto A = Matrix(4, 4, {
        2.9766, 0.3945, 0.4198, 1.1159,
        0.3945, 2.7328, -0.3097, 0.1129,
        0.4198, -0.3097, 2.5675, 0.6079,
        1.1159, 0.1129, 0.6079, 1.7231
    });

    // auto result = EIG::qr_iteration(A);
    // std::cout << result.first << "\n";
    // std::cout << "----------------\n";
    // std::cout << result.second << "\n";
    // std::cout << "----------------\n";

    // auto eigenvec = A * result.second.slice(A.get_num_rows(), 1, 0, 1);
    // std::cout << eigenvec << "\n";

    auto integral = NONLIN::midpoint_integral(0, 1, f, 0.001);
    std::cout << "Integral estimate: " << integral << "\n";   
}
