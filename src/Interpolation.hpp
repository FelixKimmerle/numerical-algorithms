#pragma once

#include <vector>

#include "Format.hpp"
#include "Polynom.hpp"

struct Point
{
    double x, y;
    Point(double x, double y) : x(x), y(y) {}
};

// Lagrange interpolation via the Lagrange polynomials (script 4.3.3)
Polynomial lagrange_interpolation(const std::vector<Point> &points);

// Aitken-Neville scheme with polynomials (script Alg. 4.12)
Polynomial aitken_neville(const std::vector<Point> &points);

// Newton's divided differences [x_0]f, [x_0,x_1]f, ..., [x_0,...,x_n]f (script Alg. 4.17)
std::vector<double> newton(const std::vector<Point> &points);

// Newton representation (script 4.44) multiplied out into the monomial basis
Polynomial get_newton_polynom(const std::vector<double> &newton_factors, const std::vector<Point> &points);

// Horner-like evaluation of the Newton representation at x (script Alg. 4.16)
double evaluate_newton(const std::vector<double> &newton_factors, const std::vector<Point> &points, double x);

// Inverse discrete Fourier transform v = M * c (script Def. 4.30)
std::vector<Complex> IDFT(const std::vector<Complex> &c);

// Discrete Fourier transform c = 1/n * conj(IDFT(conj(v))) (script 4.151)
std::vector<Complex> DFT(const std::vector<Complex> &v);
