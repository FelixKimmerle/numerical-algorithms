#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "Colors.hpp"
#include "Interpolation.hpp"
#include "Matrix.hpp"

namespace
{
std::string cell(double value)
{
    return format(value);
}

std::string cell(const Polynomial &value)
{
    std::ostringstream ss;
    ss << value;
    return ss.str();
}

// Prints the triangular scheme of Aitken-Neville / divided differences after computing entry (i, k).
// The new entry is red, the two entries it was computed from are green and the two x values used are highlighted.
// Only entries that have been computed so far are shown.
template <class T>
void print_scheme(const std::vector<std::vector<T>> &table, const std::vector<Point> &points, int i, int k)
{
    const int n = table.size();
    std::vector<size_t> widths(n, 0);
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col <= row; col++)
        {
            widths[col] = std::max(widths[col], cell(table[row][col]).size());
        }
    }

    for (int row = 0; row < n; row++)
    {
        std::ostringstream x;
        x << std::setw(8) << format(points[row].x);
        if (row == i || row == i - k)
        {
            std::cout << BLU << WHTB << x.str() << reset;
        }
        else
        {
            std::cout << x.str();
        }
        std::cout << " |";
        for (int col = 0; col <= row; col++)
        {
            const bool computed = col < k || (col == k && row <= i);
            if (!computed)
            {
                break;
            }
            std::ostringstream entry;
            entry << std::setw(widths[col]) << cell(table[row][col]);
            if (row == i && col == k)
            {
                std::cout << " " << RED << entry.str() << reset << " |";
            }
            else if ((row == i && col == k - 1) || (row == i - 1 && col == k - 1))
            {
                std::cout << " " << GRN << entry.str() << reset << " |";
            }
            else
            {
                std::cout << " " << entry.str() << " |";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Lagrange polynomial L_j: 1 in x_j and 0 in all other nodes
Polynomial lagrange(int j, const std::vector<Point> &points)
{
    Polynomial result({1});
    for (int i = 0; i < static_cast<int>(points.size()); i++)
    {
        if (i != j)
        {
            result = result * (Polynomial({-points[i].x, 1}) / (points[j].x - points[i].x));
        }
    }
    return result;
}
} // namespace

Polynomial lagrange_interpolation(const std::vector<Point> &points)
{
    Polynomial result({0});
    for (int j = 0; j < static_cast<int>(points.size()); j++)
    {
        Polynomial l = lagrange(j, points);
        std::cout << "L_" << j << "(x) = " << l << std::endl;
        result = result + l * points[j].y;
    }
    std::cout << "p(x) = sum of y_j * L_j(x)" << std::endl;
    return result;
}

Polynomial aitken_neville(const std::vector<Point> &points)
{
    const int n = points.size();
    std::vector<std::vector<Polynomial>> table(n, std::vector<Polynomial>(n, Polynomial({0})));

    for (int i = 0; i < n; i++)
    {
        table[i][0] = Polynomial({points[i].y});
    }

    for (int k = 1; k < n; k++)
    {
        for (int i = k; i < n; i++)
        {
            // P_{i,k} = (x - x_{i-k}) / (x_i - x_{i-k}) * P_{i,k-1} + (x_i - x) / (x_i - x_{i-k}) * P_{i-1,k-1}  (4.32)
            const double denominator = points[i].x - points[i - k].x;
            table[i][k] = (Polynomial({-points[i - k].x, 1}) / denominator) * table[i][k - 1] +
                          (Polynomial({points[i].x, -1}) / denominator) * table[i - 1][k - 1];
            std::cout << "P_" << i << "," << k << " = " << table[i][k] << std::endl;
            print_scheme(table, points, i, k);
        }
    }
    return table[n - 1][n - 1];
}

std::vector<double> newton(const std::vector<Point> &points)
{
    const int n = points.size();
    std::vector<std::vector<double>> table(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; i++)
    {
        table[i][0] = points[i].y;
    }

    for (int k = 1; k < n; k++)
    {
        for (int i = k; i < n; i++)
        {
            // [x_{i-k},...,x_i]f = ([x_{i-k+1},...,x_i]f - [x_{i-k},...,x_{i-1}]f) / (x_i - x_{i-k})  (4.45)
            table[i][k] = (table[i][k - 1] - table[i - 1][k - 1]) / (points[i].x - points[i - k].x);
            std::cout << "[x_" << i - k << (k == 1 ? ",x_" : ",...,x_") << i << "]f = (" << format(table[i][k - 1]) << " - "
                      << format_operand(table[i - 1][k - 1]) << ") / (" << format(points[i].x) << " - "
                      << format_operand(points[i - k].x) << ") = " << format(table[i][k]) << std::endl;
            print_scheme(table, points, i, k);
        }
    }

    // The diagonal holds [x_0]f, [x_0,x_1]f, ..., [x_0,...,x_n]f
    std::vector<double> result(n, 0.0);
    for (int i = 0; i < n; i++)
    {
        result[i] = table[i][i];
    }
    return result;
}

Polynomial get_newton_polynom(const std::vector<double> &coefficients, const std::vector<Point> &points)
{
    Polynomial result({0});
    Polynomial base({1}); // (x - x_0) * ... * (x - x_{k-1})
    std::cout << "p(x) = ";
    for (int k = 0; k < static_cast<int>(coefficients.size()); k++)
    {
        result = result + base * coefficients[k];
        if (k == 0)
        {
            std::cout << format(coefficients[k]);
        }
        else
        {
            std::cout << (coefficients[k] < 0 ? " - " : " + ") << format(std::abs(coefficients[k]));
        }
        for (int i = 0; i < k; i++)
        {
            std::cout << " * " << format_factor(points[i].x);
        }
        base = base * Polynomial({-points[k].x, 1});
    }
    std::cout << std::endl;
    return result;
}

double evaluate_newton(const std::vector<double> &coefficients, const std::vector<Point> &points, double x)
{
    const int n = coefficients.size();
    double p = coefficients[n - 1];
    std::cout << "p = " << format(p) << std::endl;
    for (int k = n - 2; k >= 0; k--)
    {
        const double next = p * (x - points[k].x) + coefficients[k];
        std::cout << "p = " << format_operand(p) << " * (" << format(x) << " - " << format_operand(points[k].x)
                  << ") + " << format_operand(coefficients[k]) << " = " << format(next) << std::endl;
        p = next;
    }
    return p;
}

std::vector<Complex> IDFT(const std::vector<Complex> &c)
{
    const int n = c.size();
    // M = (omega_n^(jk)) with omega_n = e^(2 pi i / n)  (4.143, 4.150)
    Matrix<Complex> M(n, n);
    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < n; k++)
        {
            M(j, k) = std::polar(1.0, 2.0 * M_PI * ((j * k) % n) / n);
        }
    }
    std::cout << "M =" << std::endl << M;

    Matrix<Complex> v = M * Matrix<Complex>::column(c);
    std::vector<Complex> result(n);
    for (int j = 0; j < n; j++)
    {
        result[j] = v(j, 0);
    }
    return result;
}

std::vector<Complex> DFT(const std::vector<Complex> &v)
{
    const int n = v.size();
    Matrix<Complex> conjugated = Matrix<Complex>::column(v).get_conjugated();
    std::cout << "conj(v) =" << std::endl << conjugated.transpose();

    std::vector<Complex> conjugatedValues(n);
    for (int j = 0; j < n; j++)
    {
        conjugatedValues[j] = conjugated(j, 0);
    }
    Matrix<Complex> idft = Matrix<Complex>::column(IDFT(conjugatedValues));
    std::cout << "IDFT(conj(v)) =" << std::endl << idft.transpose();

    Matrix<Complex> c = idft.get_conjugated() * Complex(1.0 / n);
    std::cout << "c = 1/" << n << " * conj(IDFT(conj(v)))" << std::endl;

    std::vector<Complex> result(n);
    for (int j = 0; j < n; j++)
    {
        result[j] = c(j, 0);
    }
    return result;
}
