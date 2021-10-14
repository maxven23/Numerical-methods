import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import *

plt.rcdefaults()

plt.style.use('classic')
plt.style.use('dark_background')

def init__func(x):
    return (math.sinh(1 / (1 + x*x)))

    
a = -3
b =  3 

n = int(input("Enter the number of divisions:  "))
h = (b - a) / n

dtype = input("Enter a differential type (RIGHT, CENTER, LEFT):  ")
#dtype = "CENTER"
title = ""
if (dtype == "CENTER"):
    title = "Центральные разности"
elif (dtype == "RIGHT"):
    title = "Правые разности"
elif (dtype == "LEFT"):
    title = "Левые разности"


#X = np.linspace(a, b, num=n)
lenght = 0
# Differetnial function
def func__double__dif(X, dif__type, func = init__func):
    F = []
    Y = [func(i) for i in X]
    length = len(X)
    for i in range(2):
        F = []
        if dif__type == "RIGHT":
            for i in range(length - 1):
                f = (Y[i + 1] - Y[i]) / h
                F.append(f)
            F.append(F[n - 2])
        elif dif__type == "CENTER":
            F.append((Y[2] - Y[0]) / (2 * h))
            for i in range(1, length - 1):
                f = (Y[i + 1] - Y[i - 1]) / (2 * h)
                F.append(f)
            F.append(F[length - 2])
        elif dif__type == "LEFT":
            F.append((Y[1] - Y[0]) / h)
            for i in range(1, length):
                f = (Y[i] - Y[i - 1]) / h
                F.append(f)
        else:
            print("SYNTAX ERROR: Unknown differential type")
        Y = F
        
    return F

#-------------------------------------------------------

init__func__x = np.linspace(a, b, num=1000)
init__func__y = [init__func(i) for i in init__func__x]

x = Symbol('x')
func = (sinh(1 / (1 + x*x)))

func = func.diff()
func = func.diff()
f = lambdify(x, func, 'numpy') 

diff__func__y = [f(i) for i in init__func__x]

#-------------------------------------------------------
# SECOND POW
def double__diff(h, n, func=init__func):
    F = []
    h = (b - a) / n
    X = np.linspace((a - h), (b + h), num=int(n + 2))
    Y = [func(i) for i in X]
    for i in range(1, len(Y) - 1):
        f = (Y[i + 1] - 2 * Y[i] + Y[i - 1]) / h**2
        F.append(f)
    return F

def comp__func(elem, h, power):
    return abs(elem) < h**power

def optimization(power):
    nb = 13
    h = (b - a) / nb
    
    #print(h**power)
    
    Xs = np.linspace(a, b, num=nb)
    
    CURR = double__diff(h, nb, init__func)
    CALC = [f(i) for i in Xs]
    RES = [CALC[i] - CURR[i] for i in range(len(CALC))]
    
    comp = [comp__func(i, h, power) for i in RES]
    
    while(not all(comp)):
        Xs = np.linspace(a, b, num=nb)
        
        CURR = double__diff(h, nb, init__func)
        CALC = [f(i) for i in Xs]
        RES = [CALC[i] - CURR[i] for i in range(len(CALC))]
        
        #print(RES)
        #print(comp)
        #print(h**power)
        
        nb += 1
        h = (b - a) / nb
        comp = [comp__func(i, h, power) for i in RES]
    #print(RES)
    #print([comp__func(i, h, power) for i in RES])
    #print(h**power)
    
    return nb

#n = optimization(2)
#print("Minimum N: " + str(n))

X = np.linspace(a, b, num=n)



fig1 = plt.subplot(211)
plt.title(title)
fig1.plot(init__func__x, diff__func__y, linewidth = 1, color="yellow")
fig1.plot(X, double__diff(h, n), "o", X, double__diff(h, n), linewidth = 1, color="white")

plt.grid(True)

plt.subplot(212)

line2 = plt.plot(init__func__x, init__func__y, color="aqua")

plt.show()