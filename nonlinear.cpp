#include "nonlinear.h"
#include <cmath>

std::pair<float, float> NONLIN::interval_bisection(float a, float b, std::function<float(float)> f, float tol) {
    while((b - a) > tol) {
        auto m = a + (b - a) / 2;
        if((f(a) > 0 && f(m) > 0) || (f(a) < 0 && f(m) < 0)) {
            a = m;
        } else {
            b = m;
        }
    }

    return std::make_pair(a, b);
}

float NONLIN::secant_method(float a, float b, std::function<float(float)> f, float tol) {
    while(fabs(f(a)) > tol) {
        auto t = b - f(b) * ((b - a) / (f(b) - f(a)));
        b = a;
        a = t;
    }
    return b;
}

float NONLIN::midpoint_integral(float a, float b, std::function<float(float)> f, float h) {
    int k = (b - a) / h;
    float I = 0;

    for(int i = 1; i <= k; ++i) {
        float m = (h / 2) * (f(a + i * h) + f(a + (i - 1) * h));
        I += m;
    }

    return I;
}
