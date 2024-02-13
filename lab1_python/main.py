import numpy as np


def show_vector(s):
    for i in range(0, len(s)):
        print(s[i], end=" ")
    print()


def fill_vectors(A):
    n = len(A)
    a = []
    b = []
    c = []
    for i in range(0, n):
        b.append(A[i][i])
    for i in range(1, n):
        a.append(A[i][i - 1])
    for i in range(1, n):
        c.append(A[i - 1][i])
    return a, b, c


def process_algorithm(a, b, c, d):
    n = len(d)
    alpha = [0] * n
    beta = [0] * n

    alpha[0] = -c[0] / b[0]
    beta[0] = d[0] / b[0]
    for i in range(1, n - 1):
        alpha[i] = -c[i] / (alpha[i - 1] * a[i - 1] + b[i])
        beta[i] = (d[i] - a[i - 1] * beta[i - 1]) / (alpha[i - 1] * a[i - 1] + b[i])

    x = [0] * n
    # x[n-1] = beta[n-1]
    x[n-1] = (d[n-1] - a[n-2] * beta[n-2]) / (a[n-2] * alpha[n-2] + b[n-1])
    for i in range(n - 2, -1, -1):
        x[i] = alpha[i] * x[i + 1] + beta[i]
    return x





def find_res(A, x):
    n = len(A)
    d = [0] * n
    for i in range(0, n):
        for j in range(0, n):
            d[i] += x[j] * A[i][j]
    return d


def main():
    n = 4
    # A = [[4, 1, 0, 0],
    #      [1, 4, 1, 0],
    #      [0, 1, 4, 1],
    #      [0, 0, 1, 4]]
    # d = [5, 6, 6, 5]
    A = [[1, 2, 0, 0],
         [2, -1, -1, 0],
         [0, 1, -1, 1],
         [0, 0, 1, 1]]
    d = [5, -3, 3, 7]
    #for i in range(1, n + 1):
    #    A[i][1:n + 1] = list(map(float, input().split()))
    #d[1:n + 1] = list(map(float, input().split()))
    print(A)
    print("d:", d)
    a, b, c = fill_vectors(A)
    show_vector(a)
    show_vector(b)
    show_vector(c)

    x = process_algorithm(a, b, c, d)
    show_vector(x)

    d_star = find_res(A, x)
    print(d_star)


    # r = (d - d*)
    r = [0] * n
    for i in range(n):
        r[i] = d[i] - d_star[i]
    # e = A^(-1) * r
    A_rev = np.matrix(A).I
    A_rev = A_rev.tolist()
    error = find_res(A_rev, r)

    print(error)


if __name__ == "__main__":
    main()