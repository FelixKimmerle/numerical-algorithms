#include <iostream>
#include <vector>
#include <cassert>

class Matrix
{
private:
    std::vector<std::vector<double>> data;
    int rows;
    int cols;

public:
    // Constructor to create an empty matrix
    Matrix(int numRows, int numCols) : rows(numRows), cols(numCols)
    {
        data.resize(rows, std::vector<double>(cols, 0.0));
    }

    // Constructor to create a matrix from a 2D vector
    Matrix(const std::vector<std::vector<double>> &input) : data(input), rows(input.size()), cols(input[0].size()) {}

    // Getter for the number of rows
    int numRows() const
    {
        return rows;
    }

    // Getter for the number of columns
    int numCols() const
    {
        return cols;
    }

    // Overloaded () operator to set a specific element in the matrix
    double &operator()(int row, int col)
    {
        if (row >= 0 && row < rows && col >= 0 && col < cols)
        {
            return data[row][col];
        }
        throw std::out_of_range("Matrix indices out of bounds.");
    }

    // Overloaded () operator for getting a specific element from the matrix
    const double &operator()(int row, int col) const
    {
        if (row >= 0 && row < rows && col >= 0 && col < cols)
        {
            return data[row][col];
        }
        throw std::out_of_range("Matrix indices out of bounds.");
    }

    // Overloaded operator* for matrix multiplication
    Matrix operator*(const Matrix &other) const
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

    // Overloaded operator<< for printing the matrix
    friend std::ostream &operator<<(std::ostream &os, const Matrix &matrix)
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

    Matrix transpose() const
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

    static Matrix identity(int size)
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

    Matrix getUpperRightTriangle(bool includeDiagonal = true) const
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
    Matrix getLowerLeftTriangle(bool includeDiagonal = true) const
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

    Matrix getDiagonal() const
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

    Matrix operator-(const Matrix &other) const
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

    static Matrix generatePermutationMatrix(int n, int index1, int index2)
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

    std::tuple<Matrix, Matrix, Matrix> LR()
    {
        assert(rows == cols);
        int n = rows;
        Matrix L = identity(n);
        Matrix R(data);
        Matrix P = identity(n);

        for (int j = 0; j < n - 1; j++)
        {
            if (R(j, j) == 0)
            {
                for (int i = j; i < n; i++)
                {
                    if (R(i, j) != 0)
                    {
                        Matrix p = generatePermutationMatrix(n,j, i);
                        R = p * R;
                        P = p * P;
                    }
                }
            }
            for (int i = j + 1; i < n; i++)
            {
                double l = R(i, j) / R(j, j);
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
    std::tuple<Matrix, Matrix> LR_Line()
    {
        assert(rows == cols);
        int n = rows;
        Matrix L = identity(n);
        Matrix R(n, n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < i; j++)
            {
                double sum = 0;
                for (int k = 0; k < j; k++)
                {
                    sum += L(i, k) * R(k, j);
                }
                L(i, j) = ((*this)(i, j) - sum) / R(j, j);
            }
            for (int j = i; j < n; j++)
            {
                double sum = 0;
                for (int k = 0; k < i; k++)
                {
                    sum += L(i, k) * R(k, j);
                }
                R(i, j) = (*this)(i, j) - sum;
            }
        }
        return std::tuple<Matrix, Matrix>(L, R.getUpperRightTriangle());
    }
};

int main()
{
    int n = 4;
    std::vector<std::vector<double>> values = {
        {0.0, 10.0, 7.0, 4.0},
        {0.0, 12.0, 1.0, 7.0},
        {0.0, 2.0, 10.0, 12.0},
        {14.0, 4.0, 12.0, 3.0}};
    Matrix A(values);
    //srand(time(0));
    //Matrix A(4, 4);
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)
    //     {
    //         A(i, j) = static_cast<int>((rand() / (double)RAND_MAX) * 10.0);
    //     }
    // }
    {
        std::cout << "A: " << std::endl
                  << A << std::endl;
        auto res = A.LR();
        auto L = std::get<0>(res);
        auto R = std::get<1>(res);
        auto P = std::get<2>(res);
        std::cout << "L:" << std::endl
                  << L << std::endl;
        std::cout << "R:" << std::endl
                  << R << std::endl;
        std::cout << "P:" << std::endl
                  << P << std::endl;
        std::cout << "L * R:" << std::endl
                  << L * R << std::endl;
        std::cout << "P * A:" << std::endl
                  << P * A << std::endl;
        std::cout << P.transpose() * P << std::endl;
    }
    // {
    //     std::cout << "A: " << std::endl
    //               << A << std::endl;
    //     auto res = A.LR_Line();
    //     auto L = std::get<0>(res);
    //     auto R = std::get<1>(res);
    //     std::cout << "L:" << std::endl
    //               << L << std::endl;
    //     std::cout << "R:" << std::endl
    //               << R << std::endl;
    //     std::cout << "L * R:" << std::endl
    //               << L * R << std::endl;
    // }
    return 0;
}
