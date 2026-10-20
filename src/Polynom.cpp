#include "Polynom.hpp"
#include "Format.hpp"

#include <algorithm>
#include <cmath>

// Constructor to initialize a polynomial with coefficients
Polynomial::Polynomial(std::vector<double> coeffs) : coefficients(coeffs) {}

// Default constructor for the zero polynomial
Polynomial::Polynomial() : coefficients{0.0} {}

// Overload the + operator for polynomial addition
Polynomial Polynomial::operator+(const Polynomial &other) const
{
    std::vector<double> result_coeffs;

    size_t max_size = std::max(coefficients.size(), other.coefficients.size());

    for (size_t i = 0; i < max_size; ++i)
    {
        double term1 = (i < coefficients.size()) ? coefficients[i] : 0.0;
        double term2 = (i < other.coefficients.size()) ? other.coefficients[i] : 0.0;
        result_coeffs.push_back(term1 + term2);
    }

    return Polynomial(result_coeffs);
}

// Overload the * operator for polynomial multiplication
Polynomial Polynomial::operator*(const Polynomial &other) const
{
    if (coefficients.empty() || other.coefficients.empty())
    {
        return Polynomial();
    }
    std::vector<double> result_coeffs(coefficients.size() + other.coefficients.size() - 1, 0.0);

    for (size_t i = 0; i < coefficients.size(); ++i)
    {
        for (size_t j = 0; j < other.coefficients.size(); ++j)
        {
            result_coeffs[i + j] += coefficients[i] * other.coefficients[j];
        }
    }

    return Polynomial(result_coeffs);
}

// Overload the * operator for scalar multiplication
Polynomial Polynomial::operator*(double scalar) const
{
    std::vector<double> result_coeffs(coefficients.size());

    for (size_t i = 0; i < coefficients.size(); ++i)
    {
        result_coeffs[i] = coefficients[i] * scalar;
    }

    return Polynomial(result_coeffs);
}

// Overload the / operator for scalar division
Polynomial Polynomial::operator/(double scalar) const
{
    std::vector<double> result_coeffs(coefficients.size());

    for (size_t i = 0; i < coefficients.size(); ++i)
    {
        result_coeffs[i] = coefficients[i] / scalar;
    }

    return Polynomial(result_coeffs);
}

// Function to evaluate the polynomial for a given value of x
double Polynomial::evaluate(double x) const
{
    double result = 0.0;
    double x_power = 1.0;

    for (double coeff : coefficients)
    {
        result += coeff * x_power;
        x_power *= x;
    }

    return result;
}

// Function to calculate and return the derivative of the polynomial
Polynomial Polynomial::derivative() const
{
    if (coefficients.size() <= 1)
    {
        return Polynomial();
    }

    std::vector<double> derivative_coeffs;
    for (size_t i = 1; i < coefficients.size(); ++i)
    {
        derivative_coeffs.push_back(coefficients[i] * i);
    }
    return Polynomial(derivative_coeffs);
}

// Prints the polynomial with the highest power first and leaves out zero terms
std::ostream &operator<<(std::ostream &os, const Polynomial &polynomial)
{
    bool first = true;
    for (int i = static_cast<int>(polynomial.coefficients.size()) - 1; i >= 0; --i)
    {
        const double value = polynomial.coefficients[i];
        if (std::abs(value) < PRINT_EPSILON)
        {
            continue;
        }
        if (first)
        {
            os << (value < 0 ? "-" : "");
        }
        else
        {
            os << (value < 0 ? " - " : " + ");
        }
        first = false;

        // Leave out a coefficient of 1 in front of x
        if (i == 0 || std::abs(std::abs(value) - 1) >= PRINT_EPSILON)
        {
            os << format(std::abs(value));
        }
        if (i > 0)
        {
            os << "x";
        }
        if (i > 1)
        {
            os << "^" << i;
        }
    }
    if (first)
    {
        os << format(0.0);
    }
    return os;
}
