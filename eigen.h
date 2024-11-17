#ifndef EIGEN_H
#define EIGEN_H

#include "matrix.h"


namespace EIG {

// Computes the dominant eigenvalue and eigenvector through normalized power iteration
std::pair<float, Matrix> normalized_power_iteration(const Matrix& A);

// Computes the eigenvalue with the least magnitude and its associated eigenvector through inverse iteration
std::pair<float, Matrix> normalized_inverse_iteration(const Matrix& A);

};


#endif