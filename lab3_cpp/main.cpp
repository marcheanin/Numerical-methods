#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

double func(double x) {
    return std::log(x) * std::log(x) / x;
}

const double exact_value = 2.0 / 3.0;

double rect_method(double a, double b, int n) {
    double h = (b - a) / n;
    double res = 0;
    for (int i = 1; i <= n; i++) {
        res += func(a + i * h + h / 2);
    }
    return res * h;
}

double trapezoid_method(double a, double b, int n) {
    double h = (b - a) / n;
    double res = 0;
    for (int i = 1; i <= n - 1; i++) {
        res += func(a + i * h);
    }
    return ((func(a) + func(b)) / 2 + res) * h;
}

double simpson_method(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0;
    for (int i = 1; i <= n - 1; i++) {
        sum += 2 * func((a + h * i)) * ((i % 2) + 1);
    }

    return (h / 3.0) * (func(a) + func(b) + sum);
}

double richardson(double Ih, double Ih2, int p) {
    return (Ih - Ih2) / (pow(2, p) - 1);
}


int main() {
    double a = 1.0 / std::numbers::e;
    double b = std::numbers::e;
    double eps = 0.001;

    std::vector <int> vector_n;
    std::vector <double> vector_Ih2;
    std::vector <double> vector_R;

    double error = 1000;
    int n = 1;
    double Ih = 0, Ih2 = 0;
    while (std::abs(error) > eps) {
        n *= 2;
        Ih = Ih2;
        Ih2 = rect_method(a, b, n);
        error = richardson(Ih, Ih2, 2);
    }
    vector_n.push_back(n);
    vector_Ih2.push_back(Ih2);
    vector_R.push_back(error);
//    std::cout   << "Rectangle method: " << std::endl
//                << "n: " << n << std::endl
//                << "I(h/2): " <<  Ih2 << std::endl
//                << "R: " << error << std::endl
//                << "I(h/2) + R: " << Ih2 + error << std::endl
//                << "Expected value: " << exact_value << std::endl;
//
//
//    std::cout << std::endl;

    error = 1000;
    n = 1;
    Ih = 0, Ih2 = 0;
    while (std::abs(error) > eps) {
        n *= 2;
        Ih2 = Ih;
        Ih = trapezoid_method(a, b, n);
        error = richardson(Ih, Ih2, 2);
    }
    vector_n.push_back(n);
    vector_Ih2.push_back(Ih2);
    vector_R.push_back(error);
//    std::cout   << "Trapezoid method: " << std::endl
//                << "n: " << n << std::endl
//                << "I(h/2): " <<  Ih2 << std::endl
//                << "R: " << error << std::endl
//                << "I(h/2) + R: " << Ih2 + error << std::endl
//                << "Expected value: " << exact_value << std::endl;
//
//
//    std::cout << std::endl;

    error = 1000;
    n = 1;
    Ih = 0, Ih2 = 0;
    while (std::abs(error) >= eps) {
        n *= 2;
        Ih2 = Ih;
        Ih = simpson_method(a, b, n);
        error = richardson(Ih, Ih2, 4);
    }
//    std::cout   << "Simpson method: " << std::endl
//                << "n: " << n << std::endl
//                << "I(h/2): " <<  Ih2 << std::endl
//                << "R: " << error << std::endl
//                << "I(h/2) + R: " << Ih2 + error << std::endl
//                << "Expected value: " << exact_value << std::endl;
    vector_n.push_back(n);
    vector_Ih2.push_back(Ih2);
    vector_R.push_back(error);

    int width = 15;

    std::cout << '\t' <<  std::setw(width) << "Rectangle method" << '\t' << "Trapezoid method" << '\t' << "Simpson method" << std::endl;
    std::cout << "n: ";
    for (int i = 0; i < vector_n.size(); i++) {
        std::cout <<  std::setw(width) << vector_n[i];
    }
    std::cout << std::endl << "I(h/2): ";
    for (int i = 0; i < vector_Ih2.size(); i++) {
        std::cout << std::setw(width) << vector_Ih2[i];
    }
    std::cout << std::endl << "R: ";
    for (int i = 0; i < vector_R.size(); i++) {
        std::cout << std::setw(width) << vector_R[i];
    }
    std::cout << std::endl << "I(h/2) + R: ";
    for (int i = 0; i < vector_R.size(); i++) {
        std::cout << std::setw(width) << vector_R[i] + vector_Ih2[i];
    }
    std::cout << std::endl << "Expected: ";
    for (int i = 0; i < vector_n.size(); i++) {
        std::cout << std::setw(width) << exact_value;
    }

}
