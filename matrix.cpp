/*
Name: Arzaq Naufail Khan
SID: 1769326
CCID: arzaqnau
AnonID: 1000376057

Description: An implementation of the matrix class
*/
#include "matrix.h"
#include <cmath>
#include <utility>
#include <cassert>

Matrix::Matrix(size_t num_rows, size_t num_columns, float init) {
    num_of_rows = num_rows;
    num_of_columns = num_columns;

    // Allocate memory for the matrix
    matrix_data = new float*[num_of_rows];
    for(size_t i = 0; i != num_of_rows; ++i) {
        matrix_data[i] = new float[num_of_columns];
    }
    
    // Initialize to init
    for(size_t i = 0; i != num_of_rows; ++i) {
        for(size_t j = 0; j != num_of_columns; ++j) {
            matrix_data[i][j] = init;
        }
    }
}

Matrix::Matrix(size_t num_rows, size_t num_columns, float *arr_ptr) {
    num_of_rows = num_rows;
    num_of_columns = num_columns;

    matrix_data = new float*[num_of_rows];
    for(size_t i = 0; i != num_of_rows; ++i) {
        matrix_data[i] = new float[num_of_columns];
    }

    for(size_t i = 0; i != num_of_rows; ++i) {
        for(size_t j = 0; j != num_of_columns; ++j) {
            // (i * num_of_columns) determines the row
            // j determines the column
            matrix_data[i][j] = arr_ptr[(i * num_of_columns) + j];
        }
    }
}

Matrix::Matrix(size_t num_rows, size_t num_columns, std::initializer_list<float> init) : num_of_rows(num_rows), num_of_columns(num_columns) {
    matrix_data = new float*[num_of_rows];
    for(int i = 0; i != num_of_rows; ++i) {
        matrix_data[i] = new float[num_columns];
    }
    auto p = init.begin();
    for(int i = 0; i != num_of_rows; ++i) {
        for(int j = 0; j != num_of_columns; ++j) {
            if(p == init.end())
                std::cerr << "Initializer list does not match matrix dimensions\n";
            matrix_data[i][j] = *p;
            ++p;
        }
    }
}

Matrix::Matrix(const Matrix& matrix) {
    num_of_rows = matrix.get_num_rows();
    num_of_columns = matrix.get_num_columns();

    // Allocate memory for elements
    matrix_data = new float*[num_of_rows];
    for(size_t i = 0; i != num_of_rows; ++i) {
        matrix_data[i] = new float[num_of_columns];
    }
    
    // Copy over the data from the matrix
    for(size_t i = 0; i != matrix.get_num_rows(); ++i) {
        for(size_t j = 0; j != matrix.get_num_columns(); ++j) {
            matrix_data[i][j] = matrix.get_element(i, j);
        }
    }

}

Matrix& Matrix::operator=(const Matrix& matrix) {
    // Delete old matrix_data
    // Need to reallocate because the new matrix may not match the current dimensions
    for(int i = 0; i != num_of_rows; ++i)
        delete[] matrix_data[i];
    delete[] matrix_data;

    num_of_rows = matrix.get_num_rows();
    num_of_columns = matrix.get_num_columns();
    
    matrix_data = new float*[num_of_rows];
    for(int i = 0; i != num_of_rows; ++i) {
        matrix_data[i] = new float[num_of_columns];
    }

    for(int i = 0; i != num_of_rows; ++i) {
        for(int j = 0; j != num_of_columns; ++j) {
            matrix_data[i][j] = matrix[i][j];
        }
    }

    return *this;
}

Matrix::Matrix(Matrix&& matrix) {
    num_of_columns = matrix.get_num_columns();
    num_of_rows = matrix.get_num_rows();
    matrix_data = matrix.clear_matrix_data(); // Careful of order, not good design, cause of many bugs
}

Matrix& Matrix::operator=(Matrix &&matrix) {
    // Delete old matrix_data
    // Need to reallocate because the new matrix may not match the current dimensions
    for(int i = 0; i != num_of_rows; ++i)
        delete[] matrix_data[i];
    delete[] matrix_data;

    num_of_rows = matrix.get_num_rows();
    num_of_columns = matrix.get_num_columns();

    matrix_data = matrix.clear_matrix_data();

    return *this;
}

Matrix Matrix::operator+(const Matrix& rhs) const {
    Matrix result = Matrix(num_of_rows, num_of_columns, 0.0f);

    for(size_t i = 0; i != num_of_rows; ++i) {
        for(size_t j = 0; j != num_of_columns; ++j) {
            result.set(i, j, matrix_data[i][j] + rhs.get_element(i, j));
        }
    }

    return result;
}

Matrix Matrix::operator-(const Matrix& rhs) const {
    Matrix result = Matrix(num_of_rows, num_of_columns, 0.0f);

    for(size_t i = 0; i != num_of_rows; ++i) {
        for(size_t j = 0; j != num_of_columns; ++j) {
            result.set(i, j, matrix_data[i][j] - rhs.get_element(i, j));
        }
    }

    return result;
}

Matrix Matrix::operator-() const {
    Matrix result = Matrix(num_of_rows, num_of_columns, 0.0f);
    
    for(size_t i = 0; i != num_of_rows; ++i) {
        for(size_t j = 0; j != num_of_columns; ++j) {
            result.set(i, j, -matrix_data[i][j]);
        }
    }

    return result;
}

Matrix Matrix::operator*(const Matrix& rhs) const {
    // M(m, n) * M(n, p) = M(m, p) 
    Matrix result = Matrix(num_of_rows, rhs.get_num_columns(), 0.0f);

    for(size_t i = 0; i != num_of_rows; ++i) {
        for(size_t j = 0; j != rhs.get_num_columns(); ++j) {
            float res = 0.0f;
            for(size_t k = 0; k != num_of_columns; ++k) {
                res += matrix_data[i][k] * rhs.get_element(k, j);
            }
            result.set(i, j, res);
        }
    }

    return result;
}

Matrix Matrix::operator*(float scalar) const {
    Matrix result(*this);
    for(auto i = 0; i != num_of_rows; ++i) {
        for(auto j = 0; j != num_of_columns; ++j) {
            result[i][j] = result[i][j] * scalar;
        }
    }

    return result;
}

Matrix Matrix::operator/(float scalar) const {
    return *this * (1 / scalar);
}

const float* Matrix::operator[](size_t index) const {
    return matrix_data[index];
}

float* Matrix::operator[](size_t index) {
    return matrix_data[index];
}

std::ostream& operator<<(std::ostream& os, Matrix &m) {
    for(size_t i = 0; i != m.get_num_rows(); ++i) {
        for(size_t j = 0; j != m.get_num_columns(); ++j) {
            if(j != m.get_num_columns() - 1) {
                std::cout << m.get_element(i, j) << " ";
            } else {
                std::cout << m.get_element(i, j); // Deal with whitespace
            }
        }
        if(i != m.get_num_rows() - 1) {
            std::cout << std::endl;
        }
    }

    return os;
}

std::istream& operator>>(std::istream& is, Matrix &m) {
    float* i_data = new float[m.get_num_rows() * m.get_num_columns()];

    for(size_t i = 0; i != m.get_num_rows() * m.get_num_columns(); ++i) {
        std::cin >> i_data[i];
    }

    for(size_t i = 0; i != m.get_num_rows(); ++i) {
        for (size_t j = 0; j != m.get_num_columns(); ++j) {
            m.set(i, j, i_data[(i * m.get_num_columns()) + j]);
        }
    }
    
    return is;
}

Matrix Matrix::transpose() const {
    Matrix result = Matrix(num_of_columns, num_of_rows, 0.0f);
    
    for(size_t i = 0; i != num_of_columns; ++i) {
        for(size_t j = 0; j != num_of_rows; ++j) {
            result.set(i, j, matrix_data[j][i]);
        }
    }

    return result;
}

Matrix::~Matrix() {
    for(size_t i = 0; i != num_of_rows; ++i) {
        delete[] matrix_data[i];
    }
    delete[] matrix_data;
}

float Matrix::get_element(size_t m, size_t n) const {
    return matrix_data[m][n];
}

size_t Matrix::get_num_rows() const {
    return num_of_rows;
}

size_t Matrix::get_num_columns() const {
    return num_of_columns;
}

void Matrix::set(size_t m, size_t n, float val) {
    matrix_data[m][n] = val;
}

float **Matrix::clear_matrix_data() {
    auto m_data = matrix_data;
    matrix_data = nullptr;
    num_of_columns = 0;
    num_of_rows = 0;

    return m_data;
    
}


Matrix Matrix::identity(size_t n) {
        Matrix I(n, n, 0.0);
        for(auto i = 0; i != n; ++i) {
            I[i][i] = 1;
        }

        return I;
}

Matrix Matrix::slice(int n, int m, int row_start, int col_start) {
    assert(n + row_start <= num_of_rows);
    assert(m + col_start <= num_of_columns);
    auto result = Matrix(n, m, 0.0f);
    
    for(auto i = 0; i != n; ++i) {
        for(auto j = 0; j != m; ++j) {
            result[i][j] = matrix_data[row_start + i][col_start + j];
        }
    }

    return result;
}

float Matrix::euclidean_norm()
{
    auto sum = 0.0f;
    for(auto i = 0; i != num_of_rows; ++i)
        sum += pow(matrix_data[i][0], 2);
    return sqrt(sum);
}

std::pair<Matrix, Matrix> gaussian_elimination(const Matrix &B) {
    Matrix A(B);
    int n = A.get_num_columns();
    Matrix M(n, n, 0.0);
    for(auto k = 0; k != n - 1; ++k) {
        if(A[k][k] == 0) {
            return std::make_pair(Matrix(n, n, 0.0), Matrix(n, n, 0.0));
        }
        for(auto i = k + 1; i != n; ++i) {
            M[i][k] = A[i][k] / A[k][k];
        }
        for(auto j = k + 1; j != n; ++j) {
            for(auto i = k + 1; i != n; ++i) {
                A[i][j] = A[i][j] - (M[i][k] * A[k][j]);
            }
        }
    }
    
    Matrix L = Matrix::identity(n);
    for(auto j = 0; j != n; ++j) {
        for(auto i = j + 1; i != n; ++i) {
            L[i][j] = M[i][j];
        }
    }
    Matrix U(n, n, 0.0);
    for(auto j = 0; j != n; ++j) {
        for(auto i = j; i >= 0; --i) {
            U[i][j] = A[i][j];
        }
    }

    return std::make_pair(L, U);
}

Matrix cholesky_factorization(Matrix A) {
    auto n = A.get_num_rows();
    for(int k = 0; k != n; ++k) {
        A[k][k] = sqrt(A[k][k]);
        for(int i = k + 1; i != n; ++i) {
            A[i][k] = A[i][k] / A[k][k]; // Scale current column
        }
        for(int j = k + 1; j != n; ++j) {
            for(int i = j; i != n; ++i) {
                A[i][j] = A[i][j] - A[i][k] * A[j][k]; // From each remaining column subtract multiple of current column
            }
        }
    }

    for(int j = 0; j != n; ++j) {
        for(int i = j - 1; i >= 0; --i) {
            A[i][j] = 0;
        }
    }

    return A;
}

Matrix forward_substitution(const Matrix& L, const Matrix& b) {
    Matrix x(b.get_num_rows(), 1, 0.0);
    Matrix b_copy(b);
    for(auto j = 0; j != b.get_num_rows(); ++j) {
        if(L[j][j] == 0)
            break;
        x[j][0] = b_copy[j][0] / L[j][j];
        for(auto i = j + 1; i != b.get_num_rows(); ++i) {
            b_copy[i][0] = b_copy[i][0] - L[i][j] * x[j][0];
        }
    }

    return x;
}

Matrix back_substitution(const Matrix& U, const Matrix& b) {
    Matrix x(b.get_num_rows(), 1, 0.0);
    Matrix b_copy(b);
    
    for(int j = b.get_num_rows() - 1; j >= 0; --j) {
        if(U[j][j] == 0)
            break;
        x[j][0] = b_copy[j][0] / U[j][j];
        for(int i = 0; i != j; ++i) {
            b_copy[i][0] = b_copy[i][0] - U[i][j] * x[j][0];
        }
    }

    return x;
}

Matrix rank_one_update_solution(std::pair<const Matrix&, const Matrix&> LU, const Matrix& u,
const Matrix& v, const Matrix& b) {
    auto z = back_substitution(LU.second, forward_substitution(LU.first, u)); // Solves Az = u
    auto y = back_substitution(LU.second, forward_substitution(LU.first, b)); // Solves Ay = b
    
    auto scaler = (v.transpose() * y)[0][0] / (1 - (v.transpose() * z)[0][0]);
    return y + (z * scaler);
}

float matrix_1_norm(const Matrix& A) {
    int n = A.get_num_rows();
    int m = A.get_num_columns();
    float sum = 0;
    for(auto j = 0; j != m; ++j) {
        float col_sum = 0;
        for(auto i = 0; i != n; ++i) {
            col_sum += fabs(A[i][j]);
        }
        if(col_sum > sum)
            sum = col_sum;
    }
    return sum;
}

float matrix_inf_norm(const Matrix& A) {
    int n = A.get_num_rows();
    int m = A.get_num_columns();
    float sum = 0;
    for(auto i = 0; i != n; ++i) {
        float row_sum = 0;
        for(auto j = 0; j != m; ++j) {
            row_sum += fabs(A[i][j]);
        }
        if(row_sum > sum)
            sum = row_sum;
    }

    return sum;
}
