#pragma once

#include "Polynom.hpp"

struct Point
{
    double x,y;
    Point(double x, double y) : x(x), y(y) {};
};


Polynomial lagrange_interpolation(const std::vector<Point> &points);
Polynomial aitken_neville(const std::vector<Point> &points);
std::vector<double> newton(const std::vector<Point> &points);
Polynomial get_newton_polynom(const std::vector<double> &newton_factors, const std::vector<Point> &points);
double evaluate_newton(const std::vector<double> &newton_factors, const std::vector<Point> &points, double x);