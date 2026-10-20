#pragma once

#include <cmath>
#include <complex>
#include <iomanip>
#include <sstream>
#include <string>

typedef std::complex<double> Complex;

// Number of decimal places used for all printed numbers
constexpr int PRECISION = 3;

// Values this close to zero are printed as 0, so rounding noise does not show up as -0.000 or 1e-17
constexpr double PRINT_EPSILON = 1e-10;

inline std::string format(double value)
{
    if (std::abs(value) < PRINT_EPSILON)
    {
        value = 0.0;
    }
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(PRECISION) << value;
    return ss.str();
}

// For subtracting a value in printed formulas: "-1" becomes "(-1.000)" so "3 - (-1.000)" stays readable
inline std::string format_operand(double value)
{
    const std::string text = format(value);
    return text[0] == '-' ? "(" + text + ")" : text;
}

// The factor (x - a) as it would be written by hand, e.g. "(x + 1.000)" or "x"
inline std::string format_factor(double a)
{
    const std::string text = format(std::abs(a));
    if (std::abs(a) < PRINT_EPSILON)
    {
        return "x";
    }
    return a < 0 ? "(x + " + text + ")" : "(x - " + text + ")";
}

inline std::string format(const Complex &value)
{
    const double re = std::abs(value.real()) < PRINT_EPSILON ? 0.0 : value.real();
    const double im = std::abs(value.imag()) < PRINT_EPSILON ? 0.0 : value.imag();
    if (im == 0.0)
    {
        return format(re);
    }
    if (re == 0.0)
    {
        return format(im) + "i";
    }
    return format(re) + (im < 0 ? " - " : " + ") + format(std::abs(im)) + "i";
}

inline std::string format_operand(const Complex &value)
{
    const std::string text = format(value);
    return text[0] == '-' || text.find(' ') != std::string::npos ? "(" + text + ")" : text;
}
