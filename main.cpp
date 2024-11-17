#include <iostream>
#include <cmath>
#include <utility>
#include "matrix.h"
#include "linear_least_squares.h"
#include "eigen.h"
#include "nonlinear.h"

#define PI 3.14159265

float my_func(float x) {
    return pow(x, 2) - (4 * sin(x));
}


int main(void) {
    auto interval = NONLIN::interval_bisection(1, 3, my_func, 0.1);
    std::cout << interval.first << " " << interval.second << "\n";
    auto sol = NONLIN::secant_method(1, 3, my_func, 0.0001);
    std::cout << "Secant method solution: " << sol << "\n";
}
