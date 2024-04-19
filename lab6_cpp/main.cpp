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

double fi_first(double x, double y) {
    return - std::pow(df_dx(x, y), 2) - std::pow(df_dy(x,y), 2);
}

double fi_second(double x, double y) {
    return df_dxdx(x, y) * std::pow(df_dx(x, y), 2) +
                    2 * df_dxdy(x, y) * df_dx(x, y) * df_dy(x, y) +
                    df_dydy(x, y) * std::pow(df_dy(x, y), 2);
}

int main() {
    double eps = 0.001 ;
    double x_k = -1, y_k = 0; // начальные значения
    double x_check = -1.011, y_check = 0.07;
    double fi1, fi2, t, error = 100;
    int n = 0;
    while (error > eps) {
        n++;
        fi1 = fi_first(x_k, y_k);
        fi2 = fi_second(x_k, y_k);
        t = - fi1 / fi2;
        x_k -= t * df_dx(x_k, y_k);
        y_k -= t * df_dy(x_k, y_k);
        error = std::max(std::abs(df_dx(x_k, y_k)), std::abs(df_dy(x_k, y_k)));
        std::cout << x_k << " " << y_k  << " " << error << std::endl;
    }
    std::cout << "Iterations: " << n << std::endl;
    std::cout << "Result: " << f(x_k, y_k) << std::endl;
    std::cout << "Analytical: " << f(x_check, y_check) << y_check << std::endl;
    std::cout << "Diff: " << std::abs(f(x_k, y_k) - f(x_check, y_check));
    return 0;
}
