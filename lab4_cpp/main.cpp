#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <iomanip>

using namespace Eigen;
using namespace std;

//double f(double x) {
//    return exp(x);
//}

//VectorXd find_lambd2(const MatrixXd& A, const VectorXd& b) {
//    return A.inverse() * b;
//}

VectorXd find_lambd(const MatrixXd& A, const VectorXd& b) {
    int m = A.rows();
    MatrixXd T = MatrixXd::Zero(m, m);
    VectorXd x = VectorXd::Zero(m);
    VectorXd y = VectorXd::Zero(m);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < i; ++j) {
            double sum = 0.0;
            for (int k = 0; k < j; ++k) {
                sum += T(i, k) * T(j, k);
            }
            T(i, j) = (A(i, j) - sum) / T(j, j);
        }

        // прямой ход
        double sum = 0.0;
        for (int k = 0; k < i; ++k) {
            sum += pow(T(i, k), 2);
        }
        T(i, i) = sqrt(A(i, i) - sum);

        // обратный ход
        double sum_y = 0.0;
        for (int k = 0; k < i; ++k) {
            sum_y += T(i, k) * y(k);
        }
        y(i) = (b(i) - sum_y) / T(i, i);
    }

    for (int i = m - 1; i >= 0; --i) {
        double sum_x = 0.0;
        for (int k = i + 1; k < m; ++k) {
            sum_x += T(k, i) * x(k);
        }
        x(i) = (y(i) - sum_x) / T(i, i);
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
    vector <double> xs = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
    vector <double> ys = {2.61, 1.62, 1.17, 0.75, 0.30, 0.75, 1.03, 0.81, 0.57}; // данные 17 варианта


    // заполняем матрицу А и свободные коэффы b по формулам:
    MatrixXd A(m, m);
    VectorXd b(m);
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

    cout << "\nA:\n" << A << endl;
    cout << "\nb:\n" << b.transpose() << endl;
    VectorXd lambd = find_lambd(A, b);
    cout << "\nλ:\n" << lambd.transpose() << endl;

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
    cout << "\nСКО Δ: " << D << endl;

    // считаем относительную ошибку
    double d = 0;
    for (int k = 0; k <= n; ++k) {
        d += pow(ys[k], 2);
    }
    d = D / sqrt(d);
    cout << fixed << setprecision(4);
    cout << "\nотн. погрешность δ: " << d << endl;

    cout << "\n|    x    |   f(x)   |   z(x)   | |f - z|  |\n";
    cout << "|---------|----------|----------|----------|\n";
    for (int k = 0; k <= n; k++) {
        cout << "| " << setw(7) << xs[k] << " | " << setw(8) << ys[k] << " | " << setw(8) << z(xs[k]) << " | " << setw(8) << abs(ys[k] - z(xs[k])) << " |\n";
        if (k != n) {
            cout << "| " << setw(7) << xs[k] + 0.25 << " | " << setw(8) << "----" << " | " << setw(8) << z(xs[k] + 0.25) << " | " << setw(8) << "----" << " |\n";
        }
    }

    return 0;
}
