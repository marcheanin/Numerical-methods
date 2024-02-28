def find_c(n, h, y):
    a = [1] * n
    c = [1] * n
    b = [4] * n
    d = [0] * n

    c[n - 1] = 0
    a[0] = 0

    for i in range(0, n):
        d[i] = 3 * (y[i + 2] - 2 * y[i + 1] + y[i]) / (h ** 2)

    alpha = [0] * n
    beta = [0] * n

    alpha[0] = -c[0] / b[0]
    beta[0] = d[0] / b[0]
    for i in range(1, n):
        alpha[i] = -c[i] / (alpha[i] * a[i - 1] + b[i])
        beta[i] = (d[i] - a[i] * beta[i - 1]) / (alpha[i] * a[i] + b[i])

    x = [0] * n
    x[n - 1] = beta[n - 1]
    # x[n - 1] = (d[n - 1] - a[n - 2] * beta[n - 2]) / (a[n - 2] * alpha[n - 2] + b[n - 1])
    for i in range(n - 2, -1, -1):
        x[i] = alpha[i] * x[i + 1] + beta[i]
    return [0] + x + [0]


def find_a_b_d(n, c, y, h):
    a = y[:-1]
    b = [0] * n
    d = [0] * n

    for i in range(n - 1):
        b[i] = (y[i + 1] - y[i]) / h - h * (c[i + 1] + 2 * c[i]) / 3
        d[i] = (c[i + 1] - c[i]) / (3 * h)

    return a, b, d


def main():
    n = 8

    # x = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
    # y = [3.33, 2.30, 1.60, 1.27, 1.18, 0.99, 1.41, 0.80, 1.12]
    x = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
    y = [3.33, 2.30, 1.60, 1.27, 1.18, 0.99, 1.41, 0.80, 1.12]

    h = (x[n] - x[0]) / n

    c = find_c(n - 1, h, y)

    a, b, d = find_a_b_d(n, c, y, h)

    x_ext = []
    splines = []
    for i in range(n):
        x_ext.append(x[i])
        x_ext.append((x[i] + x[i + 1]) / 2)
    x_ext.append(x[n])

    print("a: ", a)
    print("b: ", b)
    print("c: ", c)
    print("d: ", d)

    for i in range(len(x_ext)):
        x_star = x[0] + 0.5 * h * i
        j = i // 2
        if j == n:
            j = n - 1
        splines.append(a[j] + b[j] * (x_star - x[j]) + c[j] * (x_star - x[j]) ** 2 + d[j] * (x_star - x[j]) ** 3)

    print("x_ext: ", x_ext)

    for i in range(len(x_ext)):
        print(x_ext[i], splines[i])


if __name__ == '__main__':
    main()
