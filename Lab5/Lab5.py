import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import *

#-------------------------------------------------------

plt.style.use('classic')
plt.style.use('dark_background')

#-------------------------------------------------------

x0 = 0
y0 = 1
x_start = 0
x_end = 1

h = 1e-3

#-------------------------------------------------------

def init_func(x, y):
    return x*x - 2*y # функция первой производной

#-------------------------------------------------------
    
def Euler(n, h, x=x0, y=y0):
    for i in range(n):
        y += h * init_func(x, y)
        x += h
    return y 

#-------------------------------------------------------

def Runge_Cutta(n, h, x=x0, y=y0):
    for i in range(n):
        k1 = init_func(x, y)
        k2 = init_func(x + h/2, y + h * k1 / 2)
        k3 = init_func(x + h/2, y + h * k2 / 2)
        k4 = init_func(x + h, y + h * k3)
        y += h * (k1 + 2*k2 + 2*k3 + k4) / 6
        x += h
    return y

#-------------------------------------------------------
    
A = 1
B = math.e
n = 100
h = (b - a) / n
D0 = A + h
D1 = h

y = [[A, D0], [0, D1]]

def p(x):
    return 1
def q(x):
    return 1
def f(x):
    3 * math.e**x

def get_c1():
    global n
    return (B - y[0][n]) / y[1][n]

def get_solv_y_i(i):
    return y[0][i] + get_c1() * y[1][i]

def div(a, b):
    return a / b

def Shooting_method(n, h, x=x0, y=y0):
    for i in range(n):
        y[0].append( div((h**2 * f(x[i]) - (1 - (h / 2) * p(x[i])) * y[0][i - 1] - (h**2 * q(x[i]) - 2) * y[0]), 1 + h / 2 * p(x[i])))
        y[1].append( div(( -(1 - h/2 * p(x[i])) * y[1][i - 1] - (h**2 * q(x[i]) - 2) * y[1][i]), 1 + h /2 * p(x[i])))
    
    plt.plot(x, [get_solv_y_i(i) for i in range(n + 1)])

#-------------------------------------------------------
    
def real_func(x):
    return (2*x*x - 2*x + 3*math.e**((-2)*x) + 1) / 4

#-------------------------------------------------------
    
X = [x0 + i * h for i in range(1000)]

YE = [Euler(i, h, x0, y0) for i in range(1000)]
plt.plot(X, YE, label="Euler method")

YRC = [Runge_Cutta(i, h, x0, y0) for i in range(1000)]
plt.plot(X, YRC, label="Runge-Cutta method")

Yr = [real_func(i) for i in X]
plt.plot(X, Yr, color="aqua", label="Real solution")

#Ye = [abs(Y[i] - Yr[i]) for i in range(1000)]
#plt.plot(X, Ye, label="Error")

plt.legend(loc="upper left")
plt.grid(True)

Shooting_method(n, h, x, y)