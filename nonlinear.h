#ifndef NONLINEAR_H
#define NONLINEAR_H

#include "matrix.h"
#include <functional>

namespace NONLIN {

std::pair<float, float> interval_bisection(float a, float b, std::function<float(float)> f, float tol);

float secant_method(float a, float b, std::function<float(float)> f, float tol);

float midpoint_integral(float a, float b, std::function<float(float)> f, float h);

};

#endif