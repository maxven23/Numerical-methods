import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import *

#-------------------------------------------------------

plt.style.use('classic')
plt.style.use('dark_background')

#-------------------------------------------------------

x0 = 0
y0 = 0
x_start = 0
x_end = 1

h = 1e-3

#-------------------------------------------------------

def init_func(x, y):
    return 3*x*x*y + x*x*math.e**(x*x) # функция первой производной

#-------------------------------------------------------
    
def Euler(n, h, x=x0, y=y0):
    for i in range(n):
        y += h * init_func(x, y)
        x += h
    return y 

#-------------------------------------------------------

def real_func(x):
    return x*x*x*math.e**(x*x) / 3

#-------------------------------------------------------
    
X = [x0 + i * h for i in range(1000)]

Y = [Euler(i, h, x0, y0) for i in range(1000)]
plt.plot(X, Y, label="Our solution")

Yr = [real_func(i) for i in X]
plt.plot(X, Yr, color="aqua", label="Real solution")

Ye = [abs(Y[i] - Yr[i]) for i in range(1000)]
plt.plot(X, Ye, label="Error")

plt.legend(loc="upper left")
plt.grid(True)