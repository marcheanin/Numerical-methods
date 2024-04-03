#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

double f(double x) {
    return 2 * x;
}

double check(double x) {
    return x + std::exp(x) * std::sin(x) - std::exp(x) * std::cos(x) + 1;
}

struct equation {
    std::vector <double> as;
    std::vector <double> bs;
    std::vector <double> cs;
    std::vector <double> result;
};

std::vector<double> run(const std::vector<double>& as, const std::vector<double>& bs, const std::vector<double>& cs, const std::vector<double>& d) {
    int n = d.size();
    std::vector<double> x(n, 0), v(n, 0), u(n, 0);
    v[0] = -cs[0] / bs[0];
    u[0] = d[0] / bs[0];
    for (int i = 1; i < n; ++i) {
        v[i] = -cs[i] / (as[i] * v[i - 1] + bs[i]);
        u[i] = (d[i] - as[i] * u[i - 1]) / (as[i] * v[i - 1] + bs[i]);
    }
    x[n - 1] = u[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = v[i] * x[i + 1] + u[i];
    }
    return x;
}

equation make_equation(double h, double p, double q, int n, double x_0, double a, double b) {
    std::vector<double> as, bs, cs, result;
    bs.push_back(h * h * q - 2);
    cs.push_back(1 + h / 2 * p);
    result.push_back(h * h * f(x_0 + h) - a * (1 - h / 2 * p));
    as.push_back(0);
    for (int i = 2; i < n - 1; ++i) {
        as.push_back(1 - h / 2 * p);
        bs.push_back(h * h * q - 2);
        cs.push_back(1 + h / 2 * p);
        result.push_back(h * h * f(x_0 + h * i));
    }
    as.push_back(1 - h / 2 * p);
    bs.push_back(h * h * q - 2);
    result.push_back(h * h * f(x_0 + (n - 1) * h) - b * (1 + h / 2 * p));
    cs.push_back(0);
    return {as, bs, cs, result};
}

int main() {
    double p = -2;
    double q = 2;
    int n = 10;
    double x_0 = 0;
    double x_n = 1;
    double h = (x_n - x_0) / n;
    double a = check(x_0);
    double b = check(x_n);

    auto [as, bs, cs, result] = make_equation(h, p, q, n, x_0, a, b);
    std::vector<double> x_arr(n + 1), y_true_arr(n + 1), y_ev_arr, delta, y_gun, delta_gun;
    for (int i = 0; i <= n; ++i) {
        x_arr[i] = x_0 + i * h;
        y_true_arr[i] = check(x_0 + i * h);
    }
    y_ev_arr.push_back(a);
    auto temp = run(as, bs, cs, result);
    y_ev_arr.insert(y_ev_arr.end(), temp.begin(), temp.end());
    y_ev_arr.push_back(b);
    double max_err = 0;
    for (int i = 0; i <= n; ++i) {
        delta.push_back(std::abs(y_true_arr[i] - y_ev_arr[i]));
        max_err = std::max(delta[i], max_err);
    }
    std::cout << "Max error:" << '\t' << max_err << std::endl;
    std::cout  << std::fixed << std::setprecision(6) << "Method's res" << '\t' << "Analytical res" << '\t' << "Delta" << std::endl;

    for (int i = 0; i < delta.size(); i++) {
        std::cout << y_ev_arr[i] << "\t" << y_true_arr[i] << "\t" << delta[i] << std::endl;
    }
    std::cout << std::endl;
    return 0;
}


