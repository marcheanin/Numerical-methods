#include <iostream>
#include <cmath>

double func(double x) {
    return x * x * x + x * x - 7 * x + 4;
}

double func_dx(double x) {
    return 3 * x * x + 2 * x - 7;
}

double func_dxdx(double x) {
    return 6 * x + 2;
}

double findRoot(double a, double b, double epsilon) {
    double c;
    while ((b - a) > epsilon) {
        c = (a + b) / 2;
        if (func(c) == 0.0) {
            return c;
        } else if (func(a) * func(c) < 0) {
            b = c;
        } else {
            a = c;
        }
    }
    return (a + b) / 2;
}

int sgn(double x) {
    if (x > 0) return 1;
    if (x == 0) return 0;
    return -1;
}

double newtonMethod(double a, double b, double eps) {
    double x_prev, x_cur;
    if (func(a) * func_dxdx(a) > 0) {
        x_prev = a;
    }
    else {
        x_prev = b;
    }
    x_cur = x_prev - func(x_prev) / func_dx(x_prev);
    while (func(x_cur) * func(x_cur + sgn(x_cur - x_prev) * eps) > 0) {
       double t = x_prev - func(x_prev) / func_dx(x_prev);
       x_prev = x_cur;
       x_cur = t;
    }

    return x_cur;
}

int main() {
    std::cout << findRoot(-8, 0, 0.01) << std::endl;
    std::cout << newtonMethod(-8, 0, 0.01) << std::endl;
    return 0;
}
