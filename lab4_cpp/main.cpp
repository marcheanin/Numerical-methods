#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <iomanip>


Eigen::VectorXd find_lambd(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
    auto m = A.rows();
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(m, m);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(m);
    Eigen::VectorXd y = Eigen::VectorXd::Zero(m);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < i; ++j) {
            double sum = 0.0;
            for (int k = 0; k < j; ++k) {
                sum += U(i, k) * U(j, k);
            }
            U(i, j) = (A(i, j) - sum) / U(j, j);
        }

        // прямой ход
        double sum = 0.0;
        for (int k = 0; k < i; ++k) {
            sum += pow(U(i, k), 2);
        }
        U(i, i) = sqrt(A(i, i) - sum);

        // обратный ход
        double sum_y = 0.0;
        for (int k = 0; k < i; ++k) {
            sum_y += U(i, k) * y(k);
        }
        y(i) = (b(i) - sum_y) / U(i, i);
    }

    for (int i = m - 1; i >= 0; --i) {
        double sum_x = 0.0;
        for (int k = i + 1; k < m; ++k) {
            sum_x += U(k, i) * x(k);
        }
        x(i) = (y(i) - sum_x) / U(i, i);
    }
    return x;
}

int main() {
//    double l = 0, r = 1;
//    int m = 4, n = 10;
//    double h = (r - l) / n;
//    vector<double> xs(n + 1), ys(n + 1);
//    for (int i = 0; i <= n; ++i) {
//        xs[i] = l + i * h;
//        ys[i] = f(xs[i]);
//    }

    int m = 4; // размерность системы
    int n = 8; // кол-во точек разбиения
    std::vector <double> xs = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
    std::vector <double> ys = {2.61, 1.62, 1.17, 0.75, 0.30, 0.75, 1.03, 0.81, 0.57}; // данные 17 варианта


    // заполняем матрицу А и свободные коэффы b по формулам:
    Eigen::MatrixXd A(m, m);
    Eigen::VectorXd b(m);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            A(i, j) = 0;
            for (int k = 0; k <= n; ++k) {
                A(i, j) += pow(xs[k], i + j);
            }
        }
        b(i) = 0;
        for (int k = 0; k <= n; ++k) {
            b(i) += ys[k] * pow(xs[k], i);
        }
    }

    std::cout << "\nA:\n" << A << std::endl;
    std::cout << "\nb:\n" << b.transpose() << std::endl;
    Eigen::VectorXd lambd = find_lambd(A, b);
    std::cout << "\nλ:\n" << lambd.transpose() << std::endl;

    // задаём функцию z(x) согласно найденным лямбдам
    auto z = [&lambd, m](double x) {
        double result = 0;
        for (int i = 0; i < m; ++i) {
            result += lambd(i) * pow(x, i);
        }
        return result;
    };

    // считаем абсолютную погрешность аппроксимации
    double D = 0;
    for (int k = 0; k <= m; ++k) {
        D += pow(ys[k] - z(xs[k]), 2);
    }
    D = sqrt(D) / sqrt(n);
    std::cout << "\nСКО: " << D << std::endl;

    // считаем относительную ошибку
    double d = 0;
    for (int k = 0; k <= n; ++k) {
        d += pow(ys[k], 2);
    }
    d = D / sqrt(d);
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\nотносительная погрешность: " << d << std::endl;

    std::cout << "\n|    x    |   f(x)   |   z(x)   | |f - z|  |\n";
    std::cout << "|---------|----------|----------|----------|\n";
    for (int k = 0; k <= n; k++) {
        std::cout << "| " << std::setw(7) << xs[k] << " | " << std::setw(8) << ys[k] << " | " << std::setw(8) << z(xs[k]) << " | " << std::setw(8) << abs(ys[k] - z(xs[k])) << " |\n";
        if (k != n) { // в серединах еще посчитал z(x)
            std::cout << "| " << std::setw(7) << xs[k] + 0.25 << " | " << std::setw(8) << "----" << " | " << std::setw(8) << z(xs[k] + 0.25) << " | " << std::setw(8) << "----" << " |\n";
        }
    }

    return 0;
}
