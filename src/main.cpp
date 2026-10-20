#include <iostream>
#include <string>
#include <vector>

#include "Colors.hpp"
#include "Interpolation.hpp"
#include "Matrix.hpp"
#include "Polynom.hpp"

// TODO linear least squares (script 3.5)
// TODO Hermite interpolation (script 4.3.6)
// TODO fast IDFT / FFT (script 4.5.2)

namespace
{
void heading(const std::string &title)
{
    std::cout << std::endl << BHYEL << "======== " << title << " ========" << reset << std::endl << std::endl;
}

void check_interpolation(const Polynomial &p, const std::vector<Point> &points)
{
    std::cout << "Check: ";
    for (const auto &point : points)
    {
        std::cout << "p(" << format(point.x) << ") = " << format(p.evaluate(point.x)) << "  ";
    }
    std::cout << std::endl;
}

void print_vector(const std::string &name, const std::vector<Complex> &values)
{
    std::cout << name << " = (";
    for (size_t i = 0; i < values.size(); i++)
    {
        std::cout << (i > 0 ? ", " : "") << format(values[i]);
    }
    std::cout << ")" << std::endl;
}
} // namespace

int main()
{
    // ---------------------------------------------------------------- linear systems
    {
        heading("LR decomposition with column pivoting (Alg. 3.10)");
        Matrix<double> A({
            {0, 1, 2},
            {1, 2, 1},
            {2, 6, 1},
        });
        std::cout << "A =" << std::endl << A << std::endl;
        auto [L, R, P] = A.LR();
        std::cout << "L =" << std::endl << L << std::endl;
        std::cout << "R =" << std::endl << R << std::endl;
        std::cout << "P =" << std::endl << P << std::endl;
        std::cout << "Check: P * A - L * R =" << std::endl << P * A - L * R;
    }
    {
        heading("LR decomposition row by row (Alg. 3.8)");
        Matrix<double> A({
            {2, 1, 1},
            {4, 3, 3},
            {8, 7, 9},
        });
        std::cout << "A =" << std::endl << A << std::endl;
        auto [L, R] = A.LR_Line();
        std::cout << std::endl << "L =" << std::endl << L << std::endl;
        std::cout << "R =" << std::endl << R << std::endl;
        std::cout << "Check: A - L * R =" << std::endl << A - L * R;
    }
    {
        heading("QR decomposition with Givens rotations (Alg. 3.21)");
        Matrix<double> A({
            {3, -6},
            {4, -8},
            {0, 1},
        });
        std::cout << "A =" << std::endl << A << std::endl;
        auto [Q, R] = A.QR();
        std::cout << "Q =" << std::endl << Q << std::endl;
        std::cout << "R =" << std::endl << R << std::endl;
        std::cout << "Check: A - Q * R =" << std::endl << A - Q * R << std::endl;
        std::cout << "Check: Q^T * Q =" << std::endl << Q.transpose() * Q;
    }

    // ---------------------------------------------------------------- polynomial interpolation
    // Sampled from f(x) = x^3 - 2x + 1
    const std::vector<Point> points{{0, 1}, {1, 0}, {2, 5}, {3, 22}};
    std::cout << std::endl << "Interpolation points: ";
    for (const auto &point : points)
    {
        std::cout << "(" << format(point.x) << ", " << format(point.y) << ") ";
    }
    std::cout << std::endl;
    {
        heading("Lagrange interpolation (4.3.3)");
        Polynomial p = lagrange_interpolation(points);
        std::cout << "p(x) = " << p << std::endl;
        check_interpolation(p, points);
    }
    {
        heading("Aitken-Neville scheme (Alg. 4.12)");
        Polynomial p = aitken_neville(points);
        std::cout << "p(x) = P_" << points.size() - 1 << "," << points.size() - 1 << " = " << p << std::endl;
        check_interpolation(p, points);
    }
    {
        heading("Newton's divided differences (Alg. 4.17)");
        std::vector<double> coefficients = newton(points);
        std::cout << "Divided differences [x_0..x_k]f: ";
        for (double c : coefficients)
        {
            std::cout << format(c) << " ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "Newton representation (4.44):" << std::endl;
        Polynomial p = get_newton_polynom(coefficients, points);
        std::cout << "     = " << p << std::endl;
        check_interpolation(p, points);

        heading("Evaluation of the Newton representation at x = 4 (Alg. 4.16)");
        double value = evaluate_newton(coefficients, points, 4);
        std::cout << "p(4) = " << format(value) << "   (f(4) = 4^3 - 2 * 4 + 1 = 57)" << std::endl;
    }

    // ---------------------------------------------------------------- trigonometric interpolation
    {
        const std::vector<Complex> v{1, 2, 0, 2};
        heading("Discrete Fourier transform (Def. 4.30, 4.151)");
        print_vector("v", v);
        std::vector<Complex> c = DFT(v);
        print_vector("c", c);

        heading("Inverse discrete Fourier transform (Def. 4.30)");
        print_vector("c", c);
        std::vector<Complex> back = IDFT(c);
        print_vector("IDFT(c)", back);
        std::cout << "Check: IDFT(DFT(v)) = v" << std::endl;
    }
    return 0;
}
