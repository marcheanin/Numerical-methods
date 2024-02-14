from decimal import Decimal, getcontext

getcontext().prec = 20

def show_vector(s):
    for i in range(0, len(s)):
        print(s[i], end=" ")
    print()


def process_algorithm(a, b, c, d):
    n = len(d)
    alpha = [Decimal(0.0)] * n
    beta = [Decimal(0.0)] * n

    alpha[0] = -c[0] / b[0]
    beta[0] = d[0] / b[0]
    for i in range(1, n - 1):
        alpha[i] = -c[i] / (alpha[i - 1] * a[i - 1] + b[i])
        beta[i] = (d[i] - a[i - 1] * beta[i - 1]) / (alpha[i - 1] * a[i - 1] + b[i])

    x = [Decimal(0.0)] * n
    # x[n-1] = beta[n-1]
    x[n - 1] = (d[n - 1] - a[n - 2] * beta[n - 2]) / (a[n - 2] * alpha[n - 2] + b[n - 1])
    for i in range(n - 2, -1, -1):
        x[i] = alpha[i] * x[i + 1] + beta[i]
    return x


def check_conditions(a, b, c):
    n = len(b)
    flag = True
    for i in range(1, n - 1):
        if not (abs(b[i]) >= abs(a[i - 1]) + abs(c[i])):
            flag = False
            break

    if not (abs(b[0]) / abs(c[0]) >= 1 and abs(b[n-1]) / abs(c[n-2]) >= 1):
        flag = False

    return flag


def main():
    n = 4

    a = c = [Decimal(1.0 / 3.), Decimal(1.0 / 3.), Decimal(1.0 / 3.)]
    b = [Decimal(4.0 / 3.), Decimal(4.0 / 3.), Decimal(4.0 / 3.), Decimal(4.0 / 3.)]
    d = [Decimal(5.0 / 3.), Decimal(6.0 / 3.), Decimal(6.0 / 3.), Decimal(5.0 / 3.)]

    x_expected = [Decimal(1.0), Decimal(1.0), Decimal(1.0), Decimal(1.0)]

    if not check_conditions(a, b, c):
        print("Conditions not completed")
    else:
        print("Conditions completed")

    print()
    x = process_algorithm(a, b, c, d)

    print("x: ", end="")
    show_vector(x)

    error = [Decimal(0.0)] * n
    for i in range(n):
        error[i] = abs(x[i] - x_expected[i])

    print("error:", error)


if __name__ == "__main__":
    main()
