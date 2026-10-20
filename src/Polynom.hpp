#pragma once
#include <iostream>
#include <vector>

// Polynomial in monomial basis, coefficients[i] belongs to x^i
class Polynomial
{
private:
    std::vector<double> coefficients;

public:
    // Constructor to initialize a polynomial with coefficients
    Polynomial(std::vector<double> coeffs);

    // Default constructor for the zero polynomial
    Polynomial();

    // Overload the + operator for polynomial addition
    Polynomial operator+(const Polynomial& other) const;

    // Overload the * operator for polynomial multiplication
    Polynomial operator*(const Polynomial& other) const;

    // Overload the * operator for scalar multiplication
    Polynomial operator*(double scalar) const ;

    // Overload the / operator for scalar division
    Polynomial operator/(double scalar) const ;

    // Function to evaluate the polynomial for a given value of x
    double evaluate(double x) const;

    // Function to calculate and return the derivative of the polynomial
    Polynomial derivative() const;

    // Friend function to print the polynomial
    friend std::ostream& operator<<(std::ostream& os, const Polynomial& polynomial);
};

