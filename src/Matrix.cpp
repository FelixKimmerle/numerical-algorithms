#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "Colors.hpp"
#include "Matrix.hpp"

template <class T>
Matrix<T>::Matrix(int numRows, int numCols) : rows(numRows), cols(numCols)
{
    data.resize(rows, std::vector<T>(cols, T(0)));
}

template <class T>
Matrix<T>::Matrix(const std::vector<std::vector<T>> &input)
    : data(input), rows(input.size()), cols(input.empty() ? 0 : input[0].size())
{
    for (const auto &row : data)
    {
        if (static_cast<int>(row.size()) != cols)
        {
            throw std::invalid_argument("All rows of a matrix must have the same length.");
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::column(const std::vector<T> &values)
{
    Matrix result(values.size(), 1);
    for (size_t i = 0; i < values.size(); ++i)
    {
        result.data[i][0] = values[i];
    }
    return result;
}

template <class T>
int Matrix<T>::numRows() const
{
    return rows;
}

template <class T>
int Matrix<T>::numCols() const
{
    return cols;
}

template <class T>
T &Matrix<T>::operator()(int row, int col)
{
    if (row >= 0 && row < rows && col >= 0 && col < cols)
    {
        return data[row][col];
    }
    throw std::out_of_range("Matrix indices out of bounds.");
}

template <class T>
const T &Matrix<T>::operator()(int row, int col) const
{
    if (row >= 0 && row < rows && col >= 0 && col < cols)
    {
        return data[row][col];
    }
    throw std::out_of_range("Matrix indices out of bounds.");
}

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
Matrix<T> Matrix<T>::operator*(T scalar) const
{
    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result.data[i][j] = data[i][j] * scalar;
        }
    }
    return result;
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

template <>
Matrix<Complex> Matrix<Complex>::get_conjugated() const
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
    return *this;
}

template <class T>
Matrix<T> Matrix<T>::getUpperRightTriangle(bool includeDiagonal) const
{
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

template <class T>
Matrix<T> Matrix<T>::getLowerLeftTriangle(bool includeDiagonal) const
{
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
Matrix<T> Matrix<T>::identity(int size)
{
    Matrix identityMatrix(size, size);
    for (int i = 0; i < size; ++i)
    {
        identityMatrix.data[i][i] = T(1);
    }
    return identityMatrix;
}

template <class T>
Matrix<T> Matrix<T>::generatePermutationMatrix(int n, int index1, int index2)
{
    if (index1 < 0 || index1 >= n || index2 < 0 || index2 >= n || index1 == index2)
    {
        throw std::invalid_argument("Invalid indices for permutation matrix.");
    }

    Matrix permutationMatrix = identity(n);
    std::swap(permutationMatrix.data[index1], permutationMatrix.data[index2]);
    return permutationMatrix;
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

// All step output uses 1-based indices like the script
template <class T>
std::tuple<Matrix<T>, Matrix<T>, Matrix<T>> Matrix<T>::LR() const
{
    if (rows != cols)
    {
        throw std::invalid_argument("LR decomposition needs a square matrix.");
    }
    const int n = rows;
    Matrix L = identity(n);
    Matrix R = *this;
    Matrix P = identity(n);

    for (int j = 0; j < n - 1; j++)
    {
        std::cout << BBLU << "Step j = " << j + 1 << reset << std::endl;

        // Column pivoting: choose the element with the largest absolute value in column j (script 3.48)
        int pivot = j;
        for (int i = j + 1; i < n; i++)
        {
            if (std::abs(R(i, j)) > std::abs(R(pivot, j)))
            {
                pivot = i;
            }
        }
        if (std::abs(R(pivot, j)) == 0)
        {
            std::cout << "  Column " << j + 1 << " is already zero below the diagonal (A is singular)" << std::endl;
            continue;
        }
        if (pivot != j)
        {
            std::cout << "  Pivot: |" << format(R(pivot, j)) << "| is the largest entry in column " << j + 1
                      << " -> swap rows " << j + 1 << " and " << pivot + 1 << std::endl;
            Matrix p = generatePermutationMatrix(n, j, pivot);
            R = p * R;
            P = p * P;
            // The l_{i,k} already computed have to be swapped as well
            for (int k = 0; k < j; k++)
            {
                std::swap(L(j, k), L(pivot, k));
            }
        }
        else
        {
            std::cout << "  Pivot: " << format(R(j, j)) << " (no row swap needed)" << std::endl;
        }

        for (int i = j + 1; i < n; i++)
        {
            T l = R(i, j) / R(j, j);
            L(i, j) = l;
            std::cout << "  l_" << i + 1 << j + 1 << " = " << format(R(i, j)) << " / " << format_operand(R(j, j))
                      << " = " << format(l) << "   (row " << i + 1 << " -= " << format_operand(l) << " * row "
                      << j + 1 << ")" << std::endl;
            for (int k = j; k < n; k++)
            {
                R(i, k) -= l * R(j, k);
            }
        }
        std::cout << "  A^(" << j + 2 << ") =" << std::endl << R << std::endl;
    }
    return std::tuple<Matrix, Matrix, Matrix>(L, R, P);
}

template <class T>
std::tuple<Matrix<T>, Matrix<T>> Matrix<T>::LR_Line() const
{
    if (rows != cols)
    {
        throw std::invalid_argument("LR decomposition needs a square matrix.");
    }
    const int n = rows;
    Matrix L = identity(n);
    Matrix R(n, n);
    for (int i = 0; i < n; i++)
    {
        std::cout << BBLU << "Row i = " << i + 1 << reset << std::endl;
        for (int j = 0; j < i; j++)
        {
            T sum = 0;
            for (int k = 0; k < j; k++)
            {
                sum += L(i, k) * R(k, j);
            }
            L(i, j) = ((*this)(i, j) - sum) / R(j, j);
            std::cout << "  l_" << i + 1 << j + 1 << " = (" << format((*this)(i, j)) << " - " << format_operand(sum)
                      << ") / " << format_operand(R(j, j)) << " = " << format(L(i, j)) << std::endl;
        }
        for (int j = i; j < n; j++)
        {
            T sum = 0;
            for (int k = 0; k < i; k++)
            {
                sum += L(i, k) * R(k, j);
            }
            R(i, j) = (*this)(i, j) - sum;
            std::cout << "  r_" << i + 1 << j + 1 << " = " << format((*this)(i, j)) << " - " << format_operand(sum)
                      << " = " << format(R(i, j)) << std::endl;
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
            if (std::abs(R(i, j)) == 0)
            {
                std::cout << BBLU << "a_" << i + 1 << j + 1 << " is already 0" << reset << std::endl << std::endl;
                continue;
            }
            // Eliminate a_ij by rotating rows j and i (script Alg. 3.21, r with negative sign)
            const T r = -std::sqrt(R(j, j) * R(j, j) + R(i, j) * R(i, j));
            const T c = R(j, j) / r;
            const T s = R(i, j) / r;
            std::cout << BBLU << "Eliminate a_" << i + 1 << j + 1 << " with G_" << j + 1 << i + 1 << reset
                      << ": r = " << format(r) << ", c = " << format(c) << ", s = " << format(s) << std::endl;
            Matrix rotation = givensRotation(rows, j, i, c, s);
            std::cout << "G_" << j + 1 << i + 1 << " =" << std::endl << rotation;
            R = rotation * R;
            std::cout << "A after the rotation =" << std::endl << R << std::endl;
            G = rotation * G;
        }
    }

    // The product of all rotations is Q^T
    Matrix Q = G.transpose();
    return std::tuple<Matrix, Matrix>(Q, R);
}

template <class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &matrix)
{
    // Format everything first so the columns can be aligned
    std::vector<std::vector<std::string>> cells(matrix.rows, std::vector<std::string>(matrix.cols));
    std::vector<size_t> widths(matrix.cols, 0);
    for (int i = 0; i < matrix.rows; ++i)
    {
        for (int j = 0; j < matrix.cols; ++j)
        {
            cells[i][j] = format(matrix.data[i][j]);
            widths[j] = std::max(widths[j], cells[i][j].size());
        }
    }
    for (int i = 0; i < matrix.rows; ++i)
    {
        os << "  [";
        for (int j = 0; j < matrix.cols; ++j)
        {
            os << " " << std::setw(widths[j]) << cells[i][j];
        }
        os << " ]" << std::endl;
    }
    return os;
}

template class Matrix<double>;
template class Matrix<Complex>;
template std::ostream &operator<< <double>(std::ostream &os, const Matrix<double> &matrix);
template std::ostream &operator<< <Complex>(std::ostream &os, const Matrix<Complex> &matrix);
