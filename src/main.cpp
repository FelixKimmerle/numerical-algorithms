#include <iostream>
#include <vector>
#include <cassert>
#include <iomanip>

#include "Matrix.hpp"
#include "Polynom.hpp"
#include "Interpolation.hpp"

int main()
{
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    int n = 4;
    std::vector<std::vector<double>> values = {
        {0.0, 10.0, 7.0, 4.0},
        {0.0, 12.0, 1.0, 7.0},
        {0.0, 2.0, 10.0, 12.0},
        {14.0, 4.0, 12.0, 3.0}};
    Matrix A(values);
    // srand(time(0));
    // Matrix A(4, 4);
    //  for (int i = 0; i < n; i++)
    //  {
    //      for (int j = 0; j < n; j++)
    //      {
    //          A(i, j) = static_cast<int>((rand() / (double)RAND_MAX) * 10.0);
    //      }
    //  }
    {
        std::cout << "A: " << std::endl
                  << A << std::endl;
        std::cout << "------- LR -------" << std::endl;
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
    {
        std::cout << "A: " << std::endl
                  << A << std::endl;
        std::cout << "------- QR -------" << std::endl;
        auto res = A.QR();
        auto Q = std::get<0>(res);
        auto R = std::get<1>(res);
        std::cout << "Q:" << std::endl
                  << Q << std::endl;
        std::cout << "R:" << std::endl
                  << R << std::endl;
        std::cout << "Q * R:" << std::endl
                  << Q * R << std::endl;
    }

    std::vector<Point> points{{0,3},{1,1},{2,5},{5,-10}};
    Polynomial result = lagrange_interpolation(points);
    std::cout << result << std::endl;
    std::cout << result.evaluate(0) << std::endl;
    std::cout << result.evaluate(1) << std::endl;
    std::cout << result.evaluate(2) << std::endl;
    std::cout << result.evaluate(5) << std::endl;

    result = aitken_neville(points);
    std::cout << result << std::endl;
    std::cout << result.evaluate(0) << std::endl;
    std::cout << result.evaluate(1) << std::endl;
    std::cout << result.evaluate(2) << std::endl;
    std::cout << result.evaluate(5) << std::endl;

    auto newton_result = newton(points);
    std::cout << "newton: ";
    for(auto x : newton_result)
    {
        std::cout << x << ", ";
    }
    std::cout << std::endl;
    result = get_newton_polynom(newton_result, points);
    std::cout << result << std::endl;
    std::cout << result.evaluate(0) << std::endl;
    std::cout << result.evaluate(1) << std::endl;
    std::cout << result.evaluate(2) << std::endl;
    std::cout << result.evaluate(5) << std::endl;

    std::cout << evaluate_newton(newton_result, points, 0) << std::endl;
    std::cout << evaluate_newton(newton_result, points, 1) << std::endl;
    std::cout << evaluate_newton(newton_result, points, 2) << std::endl;
    std::cout << evaluate_newton(newton_result, points, 5) << std::endl;
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
