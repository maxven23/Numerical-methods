import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import *

#-------------------------------------------------------

plt.style.use('classic')
plt.style.use('dark_background')

#-------------------------------------------------------

eps = [0.001, 1e-6, 1e-9]

a = 0
b = 10

#-------------------------------------------------------

def init_func(x):
    return math.e**(-(x - 3)*(x - 3)) + log(1 + x) - x / 2

#-------------------------------------------------------

X = np.linspace(a, b, 5000)
Y = [init_func(x) for x in X]
plt.plot(X, Y, linewidth=1, linestyle="-", color="lime", label="Исходная функция")
plt.axhline(y=0, xmin=X[0], xmax=X[len(X) - 1], color="white", linewidth=1, label="Ось оХ")
plt.grid(True)

#-------------------------------------------------------
    
def dichotomy(a, b, eps, func=init_func):
    count = 0
    if func(a) == 0: 
        return a, count
    elif func(b) == 0:
        return b, count
    else:
        while (b - a) > eps:
            x = (a + b) / 2
            if func(x)*func(b) <= 0:
                a = x
            else:
                b = x
            count += 1
        #print("Корень получен на", count, "шаге")
        return x, count
        
#-------------------------------------------------------

x = Symbol('x')
func = math.e**(-(x - 3)*(x - 3)) + log(1 + x) - x / 2
func = func.diff()
f = lambdify(x, func, 'numpy')  
   
def Newton_method(a, b, eps, func=init_func):
    count = 0
    if func(a) == 0: 
        return a
    elif func(b) == 0:
        return b
    else:
        x = (a + b) / 2
        while abs(func(x)) > eps:
            x = x - func(x) / f(x)
            count += 1
        #print("Корень получен на", count, "шаге")
        return x, count
    
#-------------------------------------------------------

def Hord_method(a, b, eps, func=init_func):
    count = 0
    if func(a) == 0: 
        return a
    elif func(b) == 0:
        return b
    else:
        x_new = 0
        x_prev = b
        x = (a + b) / 2
        while abs(func(x)) > eps:
            x = x - func(x) * (x - (a + 0.1)) / (func(x) - func(a + 0.1))
            count += 1
        #print("Корень получен на", count, "шаге")
        return x, count
    
#-------------------------------------------------------
print("\nМетод Дихотомии\n")
for i in range(3):
    x, k = dichotomy(a, b, eps[i])
    print("Корень с точностью", eps[i], "найден на", k, "(" , 1 + floor(log(abs(b - a) / eps[i], 2)),") шаге и равен:\t", x)
    
print("\nМетод Ньютона\n")
for i in range(3):
    x, k = Newton_method(a, b, eps[i])
    print("Корень с точностью", eps[i], "найден на", k, "шаге и равен:\t", x)   
print()
print("\nМетод Хорд\n")
for i in range(3):
    x, k = Hord_method(3, b, eps[i])
    print("Корень с точностью", eps[i], "найден на", k, "шаге и равен:\t", x)   
print()

root = "3.978256718867248775074"
print("Точный корень уравнения: \t\t\t\t", root)

#-------------------------------------------------------

xs1 = []
xs2 = []
x0 = (a + b) / 2
for i in range(100):
    x0 = x0 - init_func(x0) / f(x0)
    xs2.append(x0)
    
for i in range(100):
    x0 = (a + b) / 2
    if init_func(x0)*init_func(b) <= 0:
        a = x0
    else:
        b = x0
    xs1.append(x0)

ys = [0 for i in range(100)]
plt.plot(xs1, ys, "o", color="aqua", label="Метод дихотомии")
plt.plot(xs2, ys, "o", color="gold", label="Метод Ньютона")
plt.legend(loc="lower left")
plt.show()

#-------------------------------------------------------
    
a = 0
b = 10

def CoverageRate(a, b, n, func=init_func):
    x = 3.978256718867248775074
    x0 = (a + b) / 2
    for i in range(n - 2):
        x0 = x0 - func(x0) / f(x0)

    x1 = x0 - func(x0) / f(x0)
    x2 = x1 - func(x1) / f(x1)
    x3 = x2 - func(x2) / f(x2)
    
    #print(x3, x2, x1, x0)
    #print(x3 - x2, x2 - x1, x1 - x0)
    res = math.log(abs((x3 - x2) / (x2 - x1))) / math.log(abs((x2 - x1) / (x1 - x0)))
    print("Примерная скорость сходимости (метод Ньютона):", round(res, 6))
    
    x0 = (a + b) / 2
    #for i in range(n - 2):
    
    return
    
#-------------------------------------------------------
    
a = 0
b = 10

def CoverageRateHord(a, b, n, func=init_func):
    x = 3.978256718867248775074
    x0 = (a + b) / 2
    for i in range(n - 2):
        x0 = x0 - func(x0) * (x0 - (a + 0.1)) / (func(x0) - func(a + 0.1))

    x1 = x0 - func(x0) * (x0 - (a + 0.1)) / (func(x0) - func(a + 0.1))
    x2 = x1 - func(x1) * (x1 - (a + 0.1)) / (func(x1) - func(a + 0.1))
    x3 = x2 - func(x2) * (x2 - (a + 0.1)) / (func(x2) - func(a + 0.1))
    
    #print(x3, x2, x1, x0)
    #print(x3 - x2, x2 - x1, x1 - x0)
    res = math.log(abs((x3 - x2) / (x2 - x1))) / math.log(abs((x2 - x1) / (x1 - x0)))
    print("Примерная скорость сходимости (метод Хорд):", res)
    
    x0 = (a + b) / 2
    #for i in range(n - 2):
    
    return
    
#-------------------------------------------------------
print()           
CoverageRate(a, b, 5)

CoverageRateHord(3, b, 9)
