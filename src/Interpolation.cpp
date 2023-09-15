#include <iomanip>
#include <complex>

#include "Interpolation.hpp"
#include "Matrix.hpp"
#include "Colors.hpp"
Polynomial lagrange(int j, const std::vector<Point> &points)
{
    Polynomial result({1});
    for (int i = 0; i < points.size(); i++)
    {
        if (i != j)
        {
            Polynomial n = (Polynomial({-points[i].x, 1}) / (points[j].x - points[i].x));
            result = result * n;
        }
    }
    return result;
}

Polynomial lagrange_interpolation(const std::vector<Point> &points)
{
    Polynomial result({0});

    for (int j = 0; j < points.size(); j++)
    {
        Polynomial l = lagrange(j, points);
        result = result + l * points[j].y;
    }

    return result;
}

typedef std::vector<std::vector<Polynomial>> PolTable;
typedef std::vector<std::vector<double>> NumTable;

template <class T>
void print(const std::vector<std::vector<T>> &table, const std::vector<Point> &points, Point current)
{
    for (int i = 0; i < table.size(); i++)
    {
        if (current.x == i || current.x - current.y == i)
        {
            std::cout << BLU << WHTB;
        }
        std::cout << points[i].x << reset << "|";
        for (int j = 0; j < table[i].size() && j <= i; j++)
        {
            if (current.x == i && current.y == j)
            {
                std::cout << RED;
            }
            else if (current.x == i && current.y == j + 1)
            {
                std::cout << GRN;
            }
            else if (current.x == i + 1 && current.y == j + 1)
            {
                std::cout << GRN;
            }
            int w = 9;
            if (typeid(T) == typeid(Polynomial))
            {
                w = w * (j+1);
            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << table[i][j];
            std::cout << std::setw(w) << std::setfill(' ') << ss.str() << reset << "|";
        }
        std::cout << std::endl;
    }
    std::cout << "---------" << std::endl;
}

Polynomial aitken_neville(const std::vector<Point> &points)
{
    int n = points.size();
    PolTable table(points.size(), std::vector<Polynomial>(n, Polynomial({0})));

    for (int i = 0; i < n; i++)
    {
        table[i][0] = Polynomial({points[i].y});
    }

    for (int k = 1; k < n; k++)
    {
        for (int i = k; i < n; i++)
        {
            table[i][k] = (Polynomial({-points[i - k].x, 1}) / (points[i].x - points[i - k].x)) * table[i][k - 1] + (Polynomial({points[i].x, -1}) / (points[i].x - points[i - k].x)) * table[i - 1][k - 1];
            print(table, points, Point(i, k));
        }
    }
    return table[n - 1][n - 1];
}

std::vector<double> newton(const std::vector<Point> &points)
{
    int n = points.size();
    std::vector<std::vector<double>> table(points.size(), std::vector<double>(n, 0.0));

    for (int i = 0; i < n; i++)
    {
        table[i][0] = points[i].y;
    }

    for (int k = 1; k < n; k++)
    {
        for (int i = k; i < n; i++)
        {
            // std::cout << "(" << i << "," << k << ")" << " uses: " << "(" << i << "," << k-1 << ")" << " and " << "(" << i - 1 << "," << k-1 << ")" << std::endl;
            table[i][k] = (table[i][k - 1] - table[i - 1][k - 1]) / (points[i].x - points[i - k].x);
            print(table, points, Point(i, k));
        }
    }

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
    for (int k = 0; k < coefficients.size(); k++)
    {
        Polynomial base({1});
        for (int i = 0; i < k; i++)
        {
            base = base * Polynomial({-points[i].x, 1});
        }
        result = result + base * coefficients[k];
    }
    return result;
}

double evaluate_newton(const std::vector<double> &coefficients, const std::vector<Point> &points, double x)
{
    int n = coefficients.size();
    double p = coefficients[n - 1];
    for (int k = n - 2; k >= 0; k--)
    {
        p = p * (x - points[k].x) + coefficients[k];
    }
    return p;
}

using namespace std::complex_literals;

Complex omega(double n, double p)
{
    return std::pow(std::pow(M_E, (2.0 * M_PI * 1i)/n), p);
}

Matrix<Complex> IDFT(const Matrix<Complex> &c)
{
    int n = c.numCols();
    Matrix<Complex> M(n,n);
    for(int j = 0; j < n; j++)
    {
        for(int k = 0; k < n; k++)
        {
            M(j,k) = omega(n, j * k);
        }
    }
    // auto M_ = M.transpose();
    // auto X = M * M_;
    // std::cout << X << std::endl;
    //std::cout << (M * c.transpose()) * Complex(1.0/(n)) << std::endl;
    return (M * c.transpose());
}

Matrix<Complex> DFT(const Matrix<Complex> &v)
{
    auto V_ = v.transpose().get_conjugated();
    auto idft = IDFT(V_);
    auto res = idft.get_conjugated() * (1.0/v.numRows());
    return res;
}
