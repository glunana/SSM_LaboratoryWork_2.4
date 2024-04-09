import numpy as np
from scipy.integrate import quad

def function(x):
    return 1 / (3*x**2 + 6*x + 5)

def rectangle_method(function, a, b, h):
    integral = 0
    x = a
    while x < b:
        integral += function(x) * h
        x += h
    return integral

def trapezoidal_method(function, a, b, h):
    integral = (function(a) + function(b)) / 2
    x = a + h
    while x < b:
        integral += function(x)
        x += h
    integral *= h
    return integral

def simpsons_method(function, a, b, h):
    n = int((b - a) / h)
    x = np.linspace(a, b, n+1)
    y = function(x)
    integral = h/3 * (y[0] + y[-1] + 4*np.sum(y[1:-1:2]) + 2*np.sum(y[2:-1:2]))
    return integral

def runge_romberg(h1_result, h2_result, order):
    return h2_result + (h2_result - h1_result) / ((2**order) - 1)

a = -2
b = 2
h1 = 1
h2 = 0.5

methods = {
    "Метод прямокутників": rectangle_method,
    "Метод трапецій": trapezoidal_method,
    "Метод Сімпсона": simpsons_method
}

for name, method in methods.items():
    integral_h1, _ = quad(function, a, b)
    integral_h2 = method(function, a, b, h2)
    integral_h1_approx = method(function, a, b, h1)
    integral_h2_approx = method(function, a, b, h2)

    print(f"{name}:")
    print(f"Точне значення інтегралу: {integral_h1}")
    print(f"Апроксимація з h1: {integral_h1_approx}")
    print(f"Апроксимація з h2: {integral_h2_approx}")

    for i in range(3):
        integral_h1_approx = method(function, a, b, h1 / (2**i))
        integral_h2_approx = method(function, a, b, h2 / (2**i))
        refined_integral = runge_romberg(integral_h1_approx, integral_h2_approx, i+1)
        print(f"Уточнений інтеграл (Порядок {i+1}): {refined_integral}")

    print("\n")
