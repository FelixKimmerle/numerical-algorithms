#include <cmath>
#include <complex>

#include "Matrix.hpp"

// Constructor to create an empty matrix
template <class T>
Matrix<T>::Matrix(int numRows, int numCols) : rows(numRows), cols(numCols)
{
    data.resize(rows, std::vector<T>(cols, 0.0));
}

// Constructor to create a matrix from a 2D vector
template <class T>
Matrix<T>::Matrix(const std::vector<std::vector<T>> &input) : data(input), rows(input.size()), cols(input[0].size()) {}

// Getter for the number of rows
template <class T>
int Matrix<T>::numRows() const
{
    return rows;
}

// Getter for the number of columns
template <class T>
int Matrix<T>::numCols() const
{
    return cols;
}

// Overloaded () operator to set a specific element in the matrix
template <class T>
T &Matrix<T>::operator()(int row, int col)
{
    if (row >= 0 && row < rows && col >= 0 && col < cols)
    {
        return data[row][col];
    }
    throw std::out_of_range("Matrix indices out of bounds.");
}

// Overloaded () operator for getting a specific element from the matrix
template <class T>
const T &Matrix<T>::operator()(int row, int col) const
{
    if (row >= 0 && row < rows && col >= 0 && col < cols)
    {
        return data[row][col];
    }
    throw std::out_of_range("Matrix indices out of bounds.");
}

// Overloaded operator* for matrix multiplication
template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix &other) const
{
    if (cols != other.rows)
    {
        throw std::invalid_argument("Matrix dimensions are incompatible for multiplication.");
    }

    Matrix result(rows, other.cols);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < other.cols; ++j)
        {
            for (int k = 0; k < cols; ++k)
            {
                result.data[i][j] += data[i][k] * other.data[k][j];
            }
        }
    }

    return result;
}

template <class T>
Matrix<T> Matrix<T>::operator*(T other) const
{
    Matrix result(rows, cols);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            for (int k = 0; k < cols; ++k)
            {
                result.data[i][j] = data[i][k] * other;
            }
        }
    }

    return result;
}

template <class T>
Matrix<T> Matrix<T>::transpose() const
{
    Matrix result(cols, rows);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result.data[j][i] = data[i][j];
        }
    }

    return result;
}
template <class T>
Matrix<T> Matrix<T>::identity(int size)
{
    Matrix identityMatrix(size, size);

    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            if (i == j)
            {
                identityMatrix.data[i][j] = 1.0;
            }
        }
    }

    return identityMatrix;
}

template <class T>
Matrix<T> Matrix<T>::getUpperRightTriangle(bool includeDiagonal) const
{
    if (rows != cols)
    {
        throw std::invalid_argument("Matrix is not square.");
    }

    Matrix triangleMatrix(rows, cols);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            if ((includeDiagonal && i <= j) || (!includeDiagonal && i < j))
            {
                triangleMatrix(i, j) = data[i][j];
            }
        }
    }

    return triangleMatrix;
}

// Function to get the lower left triangle of the matrix
template <class T>
Matrix<T> Matrix<T>::getLowerLeftTriangle(bool includeDiagonal) const
{
    if (rows != cols)
    {
        throw std::invalid_argument("Matrix is not square.");
    }

    Matrix triangleMatrix(rows, cols);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            if ((includeDiagonal && i >= j) || (!includeDiagonal && i > j))
            {
                triangleMatrix(i, j) = data[i][j];
            }
        }
    }

    return triangleMatrix;
}

template <class T>
Matrix<T> Matrix<T>::getDiagonal() const
{
    if (rows != cols)
    {
        throw std::invalid_argument("Matrix is not square; diagonal cannot be extracted.");
    }

    Matrix diagonalMatrix(rows, 1);

    for (int i = 0; i < rows; ++i)
    {
        diagonalMatrix(i, 0) = data[i][i];
    }

    return diagonalMatrix;
}

template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix &other) const
{
    if (rows != other.rows || cols != other.cols)
    {
        throw std::invalid_argument("Matrix dimensions are not compatible for subtraction.");
    }

    Matrix result(rows, cols);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result.data[i][j] = data[i][j] - other.data[i][j];
        }
    }

    return result;
}

template <>
Matrix<std::complex<double>> Matrix<std::complex<double>>::get_conjugated() const
{
    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result.data[i][j] = std::conj(data[i][j]);
        }
    }
    return result;
}

template <>
Matrix<double> Matrix<double>::get_conjugated() const
{
    return Matrix(data);
}

template <class T>
Matrix<T> Matrix<T>::generatePermutationMatrix(int n, int index1, int index2)
{
    if (index1 < 0 || index1 >= n || index2 < 0 || index2 >= n || index1 == index2)
    {
        throw std::invalid_argument("Invalid indices for permutation matrix.");
    }

    Matrix permutationMatrix(n, n);
    permutationMatrix = identity(n); // Initialize with an identity matrix

    // Swap rows index1 and index2 to create the permutation
    for (int col = 0; col < n; ++col)
    {
        std::swap(permutationMatrix(index1, col), permutationMatrix(index2, col));
    }

    return permutationMatrix;
}

template <class T>
std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> Matrix<T>::LR() const
{
    assert(rows == cols);
    int n = rows;
    Matrix L = identity(n);
    Matrix R(data);
    Matrix P = identity(n);

    for (int j = 0; j < n - 1; j++)
    {
        if (std::abs(R(j, j)) == 0)
        {
            for (int i = j; i < n; i++)
            {
                if (std::abs(R(i, j)) != 0)
                {
                    Matrix p = generatePermutationMatrix(n, j, i);
                    R = p * R;
                    P = p * P;
                }
            }
        }
        for (int i = j + 1; i < n; i++)
        {
            T l = R(i, j) / R(j, j);
            L(i, j) = l;
            for (int k = j; k < n; k++)
            {
                R(i, k) -= l * R(j, k);
            }
            // std::cout << "R: " << std::endl << R;
        }
    }
    return std::tuple<Matrix, Matrix, Matrix>(L, R, P);
}

template <class T>
std::tuple<Matrix<T>, Matrix<T>> Matrix<T>::LR_Line() const
{
    assert(rows == cols);
    int n = rows;
    Matrix L = identity(n);
    Matrix R(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {
            T sum = 0;
            for (int k = 0; k < j; k++)
            {
                sum += L(i, k) * R(k, j);
            }
            L(i, j) = ((*this)(i, j) - sum) / R(j, j);
        }
        for (int j = i; j < n; j++)
        {
            T sum = 0;
            for (int k = 0; k < i; k++)
            {
                sum += L(i, k) * R(k, j);
            }
            R(i, j) = (*this)(i, j) - sum;
        }
    }
    return std::tuple<Matrix, Matrix>(L, R.getUpperRightTriangle());
}

template <class T>
std::tuple<Matrix<T>, Matrix<T>> Matrix<T>::QR() const
{
    Matrix R = *this;
    Matrix G = identity(rows);
    for (int j = 0; j < cols; j++)
    {
        for (int i = j + 1; i < rows; i++)
        {
            std::cout << "begin: i: " << i << " j: " << j << std::endl;
            if (std::abs(R(i, j)) != 0)
            {
                const T r = -std::sqrt(std::pow(R(j, j), 2) + std::pow(R(i, j), 2));
                const T c = R(j, j) / r;
                const T s = R(i, j) / r;
                Matrix rotation = givensRotation(rows, j, i, c, s);
                std::cout << "Rotation: \n"
                          << rotation << std::endl;
                R = rotation * R;
                std::cout << "R after rotation: \n"
                          << R << std::endl;
                G = rotation * G;
            }
        }
    }

    Matrix Q = G.transpose();
    return std::tuple<Matrix, Matrix>(Q, R);
}

template <class T>
Matrix<T> Matrix<T>::givensRotation(int n, int i, int j, T c, T s)
{
    Matrix result = identity(n);
    result(i, i) = c;
    result(i, j) = s;
    result(j, i) = -s;
    result(j, j) = c;

    return result;
}
template class Matrix<double>;
template class Matrix<std::complex<double>>;

// Overloaded operator<< for printing the matrix
template<class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &matrix)
{
    for (int i = 0; i < matrix.rows; ++i)
    {
        for (int j = 0; j < matrix.cols; ++j)
        {
            os << matrix.data[i][j] << " ";
        }
        os << std::endl;
    }
    return os;
}
