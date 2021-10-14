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
    
def dichotomy(a, b, eps, func=init_func):
    if func(a) == 0: 
        return a
    elif func(b) == 0:
        return b
    else:
        while (b - a) > eps:
            x = (a + b) / 2
            if func(x)*func(b) <= 0:
                a = x
            else:
                b = x
        return x
        
#-------------------------------------------------------

x = Symbol('x')
func = math.e**(-(x - 3)*(x - 3)) + log(1 + x) - x / 2
func = func.diff()
f = lambdify(x, func, 'numpy')  
   
def Newton_method(a, b, eps, func=init_func):
    if func(a) == 0: 
        return a
    elif func(b) == 0:
        return b
    else:
        x = (a + b) / 2
        while abs(func(x)) > eps:
            x = x - func(x) / f(x)
        return x
#-------------------------------------------------------
    
for i in range(3):
    print("Корень с точностью ", eps[i], " методом Дихотомии равен:\t", dichotomy(a, b, eps[i]))
    print("Корень с точностью ", eps[i], " методом Ньютона равен:\t", Newton_method(a, b, eps[i]))
    
print()
root = "3.978256718867248775074"
print("Точный корень уравнения: \t\t\t\t", root)