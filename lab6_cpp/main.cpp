#include <iostream>
#include <cmath>


double f(double x, double y) {
    return x * x + 2 * y * y + 2 * x + 0.3 * std::atanh(x * y);
}

double df_dx(double x, double y) {
    return 2 * x + 2 + (3 * y) / (10 + 10 * x * x * y * y);
}

double df_dy(double x, double y) {
    return 4 * y + (3 * x) / (10 + 10 * x * x * y * y);
}

double df_dxdx(double x, double y) {
    return 2 - (60 * x * y * y * y) / ((10 + 10 * x * x * y * y) * (10 + 10 * x * x * y * y));
}

double df_dxdy(double x, double y) {
    return (30 - 30 * x * x * y * y) / ((10 + 10 * x * x * y * y) * (10 + 10 * x * x * y * y));
}

double df_dydy(double x, double y) {
    return 4 - (60 * x * x * x * y) / ((10 + 10 * x * x * y * y) * (10 + 10 * x * x * y * y));
}

int main() {
    double eps = 0.001;
    double x_k = -1, y_k = 0; // начальные значения
    double x_check = -1.011, y_check = 0.07;

    while (std::max(std::abs(df_dx(x_k, y_k)), std::abs(df_dy(x_k, y_k))) > eps) {
        
    }
}
