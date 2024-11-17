#include "eigen.h"

std::pair<float, Matrix> EIG::normalized_power_iteration(const Matrix &A) {
    auto n = A.get_num_columns();
    auto x = Matrix(n, 1, 1.3);

    Matrix y(n, 1, 0.0);
    for(auto i = 0; i != 100; ++i) {
        y = A * x;
        x = y * (1/matrix_inf_norm(y));
    }

    return std::make_pair(matrix_inf_norm(y), x);
}

std::pair<float, Matrix> EIG::normalized_inverse_iteration(const Matrix &A) {
    auto n = A.get_num_columns();
    auto x = Matrix(n, 1, 1.3);

    Matrix y(n, 1, 0.0);
    auto LU = gaussian_elimination(A);
    for(auto i = 0; i != 50; ++i) {
        y = forward_substitution(LU.first, back_substitution(LU.second, x));
        x = y * (1 / matrix_inf_norm(y));
    }

    return std::make_pair(matrix_inf_norm(y), x);
}
