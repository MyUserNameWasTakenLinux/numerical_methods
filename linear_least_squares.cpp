#include "linear_least_squares.h"
#include <cmath>


Matrix LLS::least_square_normal(const Matrix& A, const Matrix& b) {
    auto ATA = A.transpose() * A;
    auto ATb = A.transpose() * b;

    auto LU = gaussian_elimination(ATA);
    return back_substitution(LU.second, forward_substitution(LU.first, ATb));
}

Matrix LLS::householder_transform(const Matrix& v, const Matrix& u) {
    auto c = v * 2 * ((v.transpose() * u)[0][0] / (v.transpose() * v)[0][0]);
    return (u - c);
}

std::pair<std::vector<Matrix>, Matrix> LLS::householder_qr_factorization(Matrix A) {
    // m x n matrix A
    auto m = A.get_num_rows();
    auto n = A.get_num_columns();

    std::vector<Matrix> householder_transforms;

    for(auto i = 0; i != n; ++i) {
        // First calculate the householder vector
        auto v = Matrix(m, 1, 0.0f);
        for(auto j = i; j != m; ++j)
            v[j][0] = A[j][i];
        auto e_i = Matrix::identity(m).slice(m, 1, 0, i);
        auto alpha = 0.0f;
        if(A[i][i] != 0)
            alpha = (-A[i][i] / std::fabs(A[i][i])) * (v.slice(m - i, 1, i, 0).euclidean_norm());
        v = v - (e_i * alpha); // Finished calculating the householder vector
        householder_transforms.push_back(v);

        for(auto j = i; j != n; ++j) {
            auto res = householder_transform(v, A.slice(m, 1, 0, j)); // Apply householder transform to column vector
            for(auto k = 0; k != m; ++k)
                A[k][j] = res[k][0];
        }
    }

    return std::make_pair(householder_transforms, A);
}

Matrix LLS::get_q_from_transforms(std::vector<Matrix> h_transforms) {
    auto I = Matrix::identity(h_transforms[0].get_num_rows());
    
    // Iterate over columns
    for(auto i = 0; i != h_transforms.size(); ++i) {
        Matrix result = I.slice(h_transforms[0].get_num_rows(), 1, 0, i);
        for(auto v = h_transforms.rbegin(); v != h_transforms.rend(); ++v) {
            result = householder_transform(*v, result);
        }
        for(auto j = 0; j != h_transforms.size(); ++j) {
            I[j][i] = result[j][0];
        }
    }

    return I.slice(h_transforms.size(), h_transforms.size(), 0, 0);
}
