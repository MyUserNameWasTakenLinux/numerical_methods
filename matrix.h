/*
Name: Arzaq Naufail Khan
SID: 1769326
CCID: arzaqnau
AnonID: 1000376057

Description: An implementation of the matrix class
*/
#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>

class Matrix {
private:
    float **matrix_data;
    size_t num_of_columns, num_of_rows;
public:
    /*
    Description: Parameterized Constructor for the matrix

    Arguments:
        num_rows (size_t): Number of rows of the matrix
        num_columns(size_t): Number of columns of the matrix
        init (float): The value you wish the elements of the matrix to be initialized to

    Returns: N/A 
    */
    Matrix(size_t num_rows, size_t num_columns, float init);

    /*
    Description: Parameterized Constructor for the matrix

    Arguments:
        num_rows (size_t): Number of rows of the matrix
        num_columns(size_t): Number of columns of the matrix
        arr_ptr (float*): Initialize the matrix with the data contained in arr_ptr

    Returns: N/A 
    */
    Matrix(size_t num_rows, size_t num_columns, float *arr_ptr);

    /*
    Description: Constructor

    Arguments:
        num_rows (size_t)
        num_columns (size_t)
        init (initializer_list<float>): Matrix data
    
    Returns: N/A
    */

    Matrix(size_t num_rows, size_t num_columns, std::initializer_list<float> init);

    /*
    Description: Copy constructor

    Arguments:
        matrix (Matrix&): Matrix to be copied
    
    Returns: N/A
    */
    Matrix(const Matrix& matrix);

    /*
    Description: Copy assignment

    Arguments:
        matrix (Matrix&): Matrix to be copied
    
    Returns: N/A
    */
    Matrix& operator=(const Matrix& matrix);

    /*
    Description: Move constructor

    Arguments:
        matrix (Matrix&&): Matrix to be moved
    
    Returns: N/A
    */
    Matrix(Matrix&& matrix);

    /*
    Description: Move assignment

    Arguments:
        matrix (Matrix&): Matrix to be moved
    
    Returns: N/A
    */
    Matrix& operator=(Matrix&& matrix);

    /*
    Description: Adds two matrices using standard matrix addition

    Arguments:
        rhs (const Matrix&): The matrix to the right of the + operator
    
    Returns:
        A reference to a matrix that contains the result of the operation
    */
    Matrix operator+(const Matrix& rhs) const;

    /*
    Description: Substracts two matrices using standard matrix subtraction

    Arguments:
        rhs (const Matrix&): The matrix to the right of the - operator
    
    Returns:
        A reference to a matrix that contains the result of the operation
    */
    Matrix operator-(const Matrix& rhs) const;

    /*
    Description: Multiplies two matrices using standard matrix multiplication

    Arguments:
        rhs (const Matrix&): The matrix to the right of the * operator
    
    Returns:
        A reference to a matrix that contains the result of the operation
    */
    Matrix operator*(const Matrix& rhs) const;

    /*
    Description: Scalar multiplication of a matrix

    Arguments:
        rhs (float): scalar
    
    Returns:
        Matrix that contains the result of the operation
    */
   Matrix operator*(float scalar) const;

    /*
    Description: Scalar division of a matrix

    Arguments:
        rhs (float): scalar
    
    Returns:
        Matrix that contains the result of the operation
    */
   Matrix operator/(float scalar) const;

    /*
    Description: Negates a matrix using standard matrix negation

    Arguments: N/A

    Returns:
        A reference to a matrix that contains the result of the operation
    */
    Matrix operator-() const;

    /*
    Description: Indexes a (const) matrix to get one if its rows

    Arguments:
        index (size_t): The row number you want a pointer to
    
    Returns:
        A pointer to the row corresponding to index
    */
    const float* operator[](size_t index) const;

    /*
    Description: Indexes a matrix to get one if its rows

    Arguments:
        index (size_t): The row number you want a pointer to
    
    Returns:
        A pointer to the row corresponding to index
    */
    float* operator[](size_t index);

    /*
    Description: Tranpose the matrix using standard matrix transposition

    Arguments: N/A

    Returns:
        A reference to the matrix tranpose
    */
    Matrix transpose() const;

    /*
    Description: Matrix Destructor

    Arguments: N/A

    Returns: N/A
    */
    ~Matrix();

    /*
    Description: Outputs the value of the element at row m and column n

    Arguments:
        m (size_t): Row Number
        n (size_t): Column Number
    Returns:
        The value of the element at row m and column n
    */
    float get_element(size_t m, size_t n) const;
    
    /*
    Description: Outputs the value of the element at row m and column n

    Arguments:
        m (size_t): Row Number
        n (size_t): Column Number
    Returns:
        The value of the element at row m and column n
    */
    float** get_elements() const;

    /*
    Description: Outputs the number of rows in the matrix

    Arguments: N/A

    Returns:
        The number of rows in the matrix
    */

    size_t get_num_rows() const;

    /*
    Description: Outputs the number of columns in the matrix

    Arguments: N/A

    Returns:
        The number of columns in the matrix
    */
    size_t get_num_columns() const;

    /*
    Description: Sets the value of the element in row m and column n to val

    Arguments:
        m (size_t): Row number
        n (size_t): Column number
        val (float): Value to set to
    
    Returns: N/A
    */
    void set(size_t m, size_t n, float val);

    /*
    Description: Clears the matrix (DOES NOT FREE)

    Arguments:
    
    Returns:
        matrix data (float **)
    */
    float** clear_matrix_data();

    
    /*
    Description: Returns the n by n identity matrix

    Arguments:
        n (size_t): matrix dimension
    
    Returns: identity
    */
    static Matrix identity(size_t n);

    /*
    Description: Returns a slice of the matrix 

    Arguments:
        n (int): number of rows for the slice
        m (int): number of columns for the slice
        row_start (int): row index to start the slice
        col_start (int): column index to start the slice
    
    Returns: matrix slice
    */
    Matrix slice(int n, int m, int row_start, int col_start);

    float euclidean_norm();
};


/*
Description: Prints elements of the matrix

Arguments:
    os (std::ostream&): The output stream
    m (Matrix&): Matrix to output
    
Returns:
    The same output stream that was passed in
*/
std::ostream& operator<<(std::ostream& os, Matrix& m);


/*
Description: Reads in elements from the input stream into the matrix

Arguments:
    is (std::istream): The input stream
    m (Matrix&): Matrix to read into
    
Returns:
    the same input stream that was passed in
*/
std::istream& operator>>(std::istream& is, Matrix &m);


std::pair<Matrix, Matrix> gaussian_elimination(const Matrix &B);


Matrix cholesky_factorization(Matrix A);


Matrix forward_substitution(const Matrix& L, const Matrix& b);


Matrix back_substitution(const Matrix& U, const Matrix& b);


Matrix rank_one_update_solution(std::pair<const Matrix&, const Matrix&> LU, const Matrix& u,
const Matrix& v, const Matrix& b);


float matrix_1_norm(const Matrix& A);


float matrix_inf_norm(const Matrix& A);


#endif