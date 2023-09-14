#pragma once

#include <vector>
#include <cassert>
#include <iostream>


class Matrix
{
private:
    std::vector<std::vector<double>> data;
    int rows;
    int cols;

public:
    // Constructor to create an empty matrix
    Matrix(int numRows, int numCols);

    // Constructor to create a matrix from a 2D vector
    Matrix(const std::vector<std::vector<double>> &input);

    // Getter for the number of rows
    int numRows() const;

    // Getter for the number of columns
    int numCols() const;

    // Overloaded () operator to set a specific element in the matrix
    double &operator()(int row, int col);

    // Overloaded () operator for getting a specific element from the matrix
    const double &operator()(int row, int col) const;

    // Overloaded operator* for matrix multiplication
    Matrix operator*(const Matrix &other) const;

    // Overloaded operator<< for printing the matrix
    Matrix transpose() const;

    static Matrix identity(int size);

    Matrix getUpperRightTriangle(bool includeDiagonal = true) const;

    // Function to get the lower left triangle of the matrix
    Matrix getLowerLeftTriangle(bool includeDiagonal = true) const;

    Matrix getDiagonal() const;

    Matrix operator-(const Matrix &other) const;

    static Matrix generatePermutationMatrix(int n, int index1, int index2);

    std::tuple<Matrix, Matrix, Matrix> LR()const;
    std::tuple<Matrix, Matrix> LR_Line()const;
    std::tuple<Matrix, Matrix> QR()const;
    static Matrix givensRotation(int n, int j, int i, double c, double s);

    friend std::ostream &operator<<(std::ostream &os, const Matrix &matrix);
};