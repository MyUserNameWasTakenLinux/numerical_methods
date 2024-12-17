#ifndef LINEAR_LEAST_SQUARES_HPP
#define LINEAR_LEAST_SQUARES_HPP

#include "matrix.h"
#include <vector>

namespace LLS {

// Solve linear least squares using the normal equations
// Returns linear least squares solution to Ax = b
Matrix least_square_normal(const Matrix& A, const Matrix& b);

// Apply the householder transform to u, where v is the householder vector
Matrix householder_transform(const Matrix& v, const Matrix& u);

// Returns the householder vectors in a std::vector<> and an upper triangular matrix
std::pair<std::vector<Matrix>, Matrix> householder_qr_factorization(Matrix A);

// Reduce QR factorization Q
Matrix get_q_from_transforms(std::vector<Matrix> h_transforms);

};

#endif