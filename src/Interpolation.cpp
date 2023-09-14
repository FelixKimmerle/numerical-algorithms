#include "Interpolation.hpp"

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

Polynomial aitken_neville(const std::vector<Point> &points)
{
    int n = points.size();
    std::vector<std::vector<Polynomial>> table(points.size(), std::vector<Polynomial>(n, Polynomial({0})));

    for (int i = 0; i < n; i++)
    {
        table[i][0] = Polynomial({points[i].y});
    }

    for (int k = 1; k < n; k++)
    {
        for (int i = k; i < n; i++)
        {
            table[i][k] = (Polynomial({-points[i - k].x, 1}) / (points[i].x - points[i - k].x)) * table[i][k - 1] + (Polynomial({points[i].x, -1}) / (points[i].x - points[i - k].x)) * table[i - 1][k - 1];
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
            //std::cout << "(" << i << "," << k << ")" << " uses: " << "(" << i << "," << k-1 << ")" << " and " << "(" << i - 1 << "," << k-1 << ")" << std::endl;
            table[i][k] = (table[i][k - 1] - table[i - 1][k - 1]) / (points[i].x - points[i - k].x);
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
