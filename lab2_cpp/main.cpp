#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <iomanip>

double func(double x) {
    return std::log(x) * std::log(x) / x;
}

std::vector<double> find_c(int n, double h, const std::vector<double>& y) {
    std::vector<double> a(n, 1), c(n, 1), b(n, 4), d(n, 0);
    c[n - 1] = 0;
    a[0] = 0;
    for (int i = 0; i < n; ++i) {
        d[i] = 3 * (y[i + 2] - 2 * y[i + 1] + y[i]) / (h * h);
    }
    std::vector<double> alpha(n, 0), beta(n, 0);
    alpha[0] = -c[0] / b[0];
    beta[0] = d[0] / b[0];
    for (int i = 1; i < n; ++i) {
        alpha[i] = -c[i] / (alpha[i - 1] * a[i - 1] + b[i]);
        beta[i] = (d[i] - a[i] * beta[i - 1]) / (alpha[i - 1] * a[i - 1] + b[i]);
    }
    std::vector<double> x(n, 0);
    x[n - 1] = beta[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
    std::vector<double> result = {0};
    result.insert(result.end(), x.begin(), x.end());
    result.push_back(0);
    return result;
}

std::tuple <std::vector<double>, std::vector<double>, std::vector<double> > find_a_b_d(int n, const std::vector<double>& c, const std::vector<double>& y, double h) {
    std::vector<double> a(y.begin(), y.end() - 1), b(n, 0), d(n, 0);
    for (int i = 0; i < n - 1; ++i) {
        b[i] = (y[i + 1] - y[i]) / h - h * (c[i + 1] + 2 * c[i]) / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h);
    }
    return {a, b, d};
}

int main() {
    int n = 31;

    std::vector <double> x, y;

    auto a_dot = 1 / std::numbers::e;
    auto b_dot = std::numbers::e;

    auto step = (b_dot - a_dot) / n;

    for (int i = 0; i <= n; i++){
        x.push_back(a_dot + step * i);
        y.push_back(func(a_dot + step * i));
    }

    for (int i = 0; i < x.size(); i++){
        std::cout << i << ": " << x[i] << " " << y[i] << std::endl;
    }

    std::cout << std::endl;

    double h = (x[n] - x[0]) / n;
    std::vector<double> c = find_c(n - 1, h, y);
    auto res = find_a_b_d(n, c, y, h);
    auto a = std::get<0>(res);
    auto b = std::get<1>(res);
    auto d = std::get<2>(res);
    std::vector<double> x_ext;
    std::vector<double> splines;
    for (int i = 0; i < n; ++i) {
        x_ext.push_back(x[i]);
        x_ext.push_back((x[i] + x[i + 1]) / 2);
    }
    x_ext.push_back(x[n]);
    std::cout << "a: ";
    for (const auto& val : a) std::cout << val << " ";
    std::cout << std::endl;
    std::cout << "b: ";
    for (const auto& val : b) std::cout << val << " ";
    std::cout << std::endl;
    std::cout << "c: ";
    for (const auto& val : c) std::cout << val << " ";
    std::cout << std::endl;
    std::cout << "d: ";
    for (const auto& val : d) std::cout << val << " ";
    std::cout << std::endl;
    for (int i = 0; i < x_ext.size(); ++i) {
        double x_star = x[0] + 0.5 * h * i;
        int j = i / 2;
        if (j == n) {
            j = n - 1;
        }
        splines.push_back(a[j] + b[j] * (x_star - x[j]) + c[j] * (x_star - x[j]) * (x_star - x[j]) + d[j] * (x_star - x[j]) * (x_star - x[j]) * (x_star - x[j]));
    }
    std::cout << "x_ext: ";
    for (const auto& val : x_ext) std::cout << val << " ";
    std::cout << std::endl << std::endl;

    int j = 0;
    for (int i = 0; i < x_ext.size(); ++i) {
        std::cout << std::fixed << std::setprecision(4) << x_ext[i] << " " << splines[i];
        if (i % 2 == 0) {
            std::cout << std::fixed << std::setprecision(4) << "\t" << y[j] << "\t" << std::abs(splines[i] - y[j]) << "\t -- node";
            j++;
        }
        else {
            std::cout << std::fixed << std::setprecision(4) << "\t" << func(x_ext[i]) << "\t" << std::abs(splines[i] - y[j]);
        }
        std::cout << std::endl;
    }
    return 0;
}


