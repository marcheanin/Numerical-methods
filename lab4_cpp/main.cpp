#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <iomanip>

using namespace Eigen;

std::vector <std::vector <int> > funcs = {{1, 3, 6}, {4, 2, 9}, {5, 8, 7}};

Eigen::VectorXd findLambd(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
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

double avg(double x1, double x2) {
    return (x1 + x2) / 2;
}

double geom(double x1, double x2) {
    return sqrt(x1 * x2);
}

double harm(double x1, double x2) {
    return 2 / (1 / x1 + 1 / x2);
}

int findK(std::vector <double> xs, std::vector <double> ys, const std::function<double(double)>& z) {
    double x_first = xs[0], x_last = xs[xs.size()-1];
    double y_first = ys[0], y_last = ys[ys.size()-1];
    std::vector <double> funcs_x = {avg(x_first, x_last),   // x_a
                                    geom(x_first, x_last),  // x_g
                                    harm(x_first, x_last)}; // x_h
    std::vector <double> funcs_y = {avg(y_first, y_last),   // y_a
                                    geom(y_first, y_last),  // y_g
                                    harm(y_first, y_last)}; // y_h
    double min = abs(z(funcs_x[0]) - funcs_y[0]);
    int min_x = 0, min_y = 0;
    for (int i = 0; i < funcs_x.size(); i++) {
        for (int j = 0; j < funcs_y.size(); j++) {
            if (std::abs(z(funcs_x[i]) - z(funcs_y[j])) < min) {
                min = std::abs(z(funcs_x[i]) - z(funcs_y[j]));
                min_x = i;
                min_y = j;
            }
        }
    }
    return funcs[min_x][min_y];
}

void fillAb (MatrixXd & A, MatrixXd & b, int m, std::vector <double> xs, std::vector <double> ys, int n) {
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
}


int main() {

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
    Eigen::VectorXd lambd = findLambd(A, b);
    std::cout << "\nλ:\n" << lambd.transpose() << std::endl;

    // задаём функцию z(x) согласно найденным лямбдам
    auto z = [&lambd, m](double x) {
        double result = 0;
        for (int i = 0; i < m; ++i) {
            result += lambd(i) * pow(x, i);
        }
        return result;
    };

    std::cout << "k: " << findK(xs, ys, z) << std::endl; // результат - 5



    Eigen::MatrixXd A_new(2, 2);
    Eigen::VectorXd b_new(2);

    std::vector <double> xs_rev (xs.size());
    for (int i = 0; i < xs.size(); i++) {
        xs_rev[i] = 1 / xs[i];
    }

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            A_new(i, j) = 0;
            for (int k = 0; k <= n; ++k) {
                A_new(i, j) += pow(xs_rev[k], i + j);
            }
        }
        b_new(i) = 0;
        for (int k = 0; k <= n; ++k) {
            b_new(i) += ys[k] * pow(xs_rev[k], i);
        }
    }

    VectorXd new_lambd = findLambd(A_new, b_new);

    std::cout << "New lambd: " <<  std::endl <<  new_lambd.transpose();
    auto z_new = [&new_lambd](double x) {
        return new_lambd[0] / x + new_lambd[1];
    };

    // считаем абсолютную погрешность аппроксимации
    double D = 0;
    for (int k = 0; k <= m; ++k) {
        D += pow(ys[k] - z_new(xs[k]), 2);
    }


    D = sqrt(D) / sqrt(n);
    std::cout << "\nСКО: " << D << std::endl;

    // считаем относительню ошибку
    double d = 0;
    for (int k = 0; k <= n; ++k) {
        d += pow(ys[k], 2);
    }
    d = D / sqrt(d);
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\nотносительная погрешность: " << d << std::endl;

    return 0;
}
