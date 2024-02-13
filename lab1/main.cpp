#include <iostream>
#include <vector>


void showVector(const std::vector <double>& s) {
    for (int i = 1; i < s.size(); i++){
        std::cout << s[i] << " ";
    }
    std::cout << std::endl;
}

void fillVectors(const std::vector <std::vector <double> >& A, std::vector <double>& a, std::vector <double>& b, std::vector <double>& c) {
    int n = A.size() - 1;
    a.push_back(0), b.push_back(0), c.push_back(0);
    for (int i = 1; i <= n; i++) {
        b.push_back(A[i][i]);
    }
    for (int i = 2; i <= n; i++) {
        a.push_back(A[i][i - 1]);
    }
    for (int i = 2; i <= n; i++){
        c.push_back(A[i - 1][i]);
    }
}

std::vector <double> processAlgorithm(const std::vector <double> &a, const std::vector <double> &b, const std::vector <double> &c, const std::vector <double> &d) {
    std::vector <double> alpha; std::vector <double> beta;
    int n = a.size();
    std::cout << n << std::endl;
    alpha.resize(n + 1), beta.resize(n + 1);

    // forward
    alpha[1] = - c[1] / b[1];
    beta[1] = d[1] / b[1];

    for (int i = 2; i <= n; i++) {
        alpha[i] = -c[i] / (alpha[i - 1] * a[i - 1] + b[i]);
        beta[i] = (d[i] - a[i - 1]*beta[i - 1]) / (alpha[i-1]*a[i-1] + b[i]);
    }

    // backward
    std::vector <double> x (n + 1);
    x[n] = beta[n];
    for (int i = n-1; i >= 1; i--){
        x[i] = alpha[i] * x[i + 1] + beta[i]    ;
    }

    return x;
}

std::vector <double> findRes(const std::vector <std::vector <double> >& A, const std::vector <double>& x) {
    int n = A.size() - 1;
    std::vector <double> d (n + 1);
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            d[i] += x[j] * A[i][j];
        }
    }
    return d;
}


int main() {
    int n;
    std::cin >> n;

    std::vector <std::vector <double> > A (n + 1, std::vector <double> (n + 1, 0));
    std::vector <double> d (n + 1);

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            std::cin >> A[i][j];
        }
    }
    for (int i = 1; i <= n; i++){
        std::cin >> d[i];
    }

    std::vector <double> a, b, c;

    fillVectors(A, a, b, c);
    showVector(a), showVector(b), showVector(c);

    std::vector <double> x = processAlgorithm(a, b, c, d);

    showVector(x);

    // count error

    auto d_star = findRes(A, x);
    std::vector <double> r;
}
/*
 4
 4 1 0 0
 1 4 1 0
 0 1 4 1
 0 0 1 4
 5 6 6 5

 4
 1 2 0 0
 2 -1 -1 0
 0 1 -1 1
 0 0 1 1

 5 -3 3 7
 */
