#pragma once

#include <iostream>
#include <tuple>
#include <vector>

#include "Format.hpp"

template <class T>
class Matrix;

template <class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &matrix);

template <class T>
class Matrix
{
private:
    std::vector<std::vector<T>> data;
    int rows;
    int cols;

public:
    // Constructor to create a zero matrix
    Matrix(int numRows, int numCols);

    // Constructor to create a matrix from a 2D vector
    Matrix(const std::vector<std::vector<T>> &input);

    // Constructor to create a column vector
    static Matrix column(const std::vector<T> &values);

    int numRows() const;
    int numCols() const;

    // Access to a specific element (0-based)
    T &operator()(int row, int col);
    const T &operator()(int row, int col) const;

    Matrix operator*(const Matrix &other) const;
    Matrix operator*(T scalar) const;
    Matrix operator-(const Matrix &other) const;

    Matrix transpose() const;
    Matrix get_conjugated() const;
    Matrix getUpperRightTriangle(bool includeDiagonal = true) const;
    Matrix getLowerLeftTriangle(bool includeDiagonal = true) const;

    // Diagonal entries as a column vector
    Matrix getDiagonal() const;

    static Matrix identity(int size);

    // Transposition matrix P_{index1,index2}: the identity with rows index1 and index2 swapped (script Def. 3.9)
    static Matrix generatePermutationMatrix(int n, int index1, int index2);

    // Givens rotation G_{i,j}(c, s) (script Def. 3.19)
    static Matrix givensRotation(int n, int i, int j, T c, T s);

    // LR decomposition with column pivoting (script Alg. 3.10), returns (L, R, P) with P * A = L * R
    std::tuple<Matrix, Matrix, Matrix> LR() const;

    // LR decomposition computed row by row, without pivoting (script Alg. 3.8), returns (L, R) with A = L * R
    std::tuple<Matrix, Matrix> LR_Line() const;

    // QR decomposition with Givens rotations (script Alg. 3.21), returns (Q, R) with A = Q * R
    std::tuple<Matrix, Matrix> QR() const;

    friend std::ostream &operator<< <>(std::ostream &os, const Matrix<T> &matrix);
};
