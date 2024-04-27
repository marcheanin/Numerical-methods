import numpy as np


def newtonMethod(func, J, x0, eps=0.01):
    x = x0
    for _ in range(100):
        y = np.linalg.solve(J(x), -func(x))
        x_new = x + y
        if np.max(np.abs(x_new - x)) < eps:
            return x_new
        x = x_new

    print("Алгоритм не сошелся за 100 итераций")
    return x


def func(x):
    return np.array([
        np.cos(x[0] - 1) + x[1] - 0.8,
        x[0] - np.cos(x[1]) - 2
    ])


def J(x):
    return np.array([
        [-np.sin(x[0] - 1), 1],
        [1, np.sin(x[1])]
    ])


x0 = np.array([3.0, 1.0])
solution = newtonMethod(func, J, x0)
print(solution)
