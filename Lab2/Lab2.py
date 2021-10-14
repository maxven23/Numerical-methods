import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import *

#-------------------------------------------------------

plt.style.use('classic')
plt.style.use('dark_background')

#-------------------------------------------------------

def init__func(x):
    return (math.sinh(1 / (1 + x*x)))

#-------------------------------------------------------
    
a = -3
b =  3 

n = int(input("Enter the number of divisions:  "))
h = (b - a) / n

#---------------Задаем тип разностей------------------

dtype = input("Enter a differential type (RIGHT, CENTER, LEFT):  ")
title = ""
if (dtype == "CENTER"):
    title = "Центральные разности"
elif (dtype == "RIGHT"):
    title = "Правые разности"
elif (dtype == "LEFT"):
    title = "Левые разности"

#-----------------Differetnial function-----------------
def func__dif(dif__type, n, func = init__func):
    F = []
    h = (b - a) / n
    if dif__type == "RIGHT":
        X = np.linspace(a, b + h, n + 1)
        Y = [func(i) for i in X]
        
        for i in range(n):
            f = (Y[i + 1] - Y[i]) / h
            F.append(f)
            
    elif dif__type == "CENTER":
        X = np.linspace(a - h, b + h, n + 2)
        Y = [func(i) for i in X]
        
        for i in range(1, n + 1):
            f = (Y[i + 1] - Y[i - 1]) / (2 * h)
            F.append(f)
            
    elif dif__type == "LEFT":
        X = np.linspace(a, b + h, n + 1)
        Y = [func(i) for i in X]
        
        for i in range(n):
            f = (Y[i] - Y[i - 1]) / h
            F.append(f)
    else:
        print("SYNTAX ERROR: Unknown differential type")
    
    return F

x = Symbol('x')
func = (sinh(1 / (1 + x*x)))

func = func.diff()
f = lambdify(x, func, 'numpy') 


#-------------------------------------------------------
"""def difference(correct, got):
    res = 0
    for i in range(len(correct)):
        res += (correct[i] - got[i])
    res /= len(correct)
    return abs(res)
"""
def difference(correct, got):
    res = [abs(correct[i] - got[i]) for i in range(len(got))]
    return max(res)

def log__graph():
    res = []
    hs = []
    for j in range(3, 1000):
        h = (b - a) / j
        X = np.linspace(a, b, num=j)
        Y1 = func__dif(dtype, j)
        Y2 = [f(i) for i in X]
        hs.append(log(h))
        res.append(log(difference(Y2, Y1)))
    return hs, res
#-------------------------------------------------------
    
init__func__x = np.linspace(a, b, num=1000)
init__func__y = [init__func(i) for i in init__func__x]
diff__func__y = [f(i) for i in init__func__x]

#-------------------------------------------------------


X = np.linspace(a, b, n)
    
fig1 = plt.subplot(211)
plt.title(title)
plt.ylabel("y'(x)")

fig1.plot(init__func__x, diff__func__y, linewidth = 1, color="yellow")
fig1.plot(X, func__dif(dtype, n), "o", X, func__dif(dtype, n), linewidth = 1, color="white")

plt.grid(True)

#-------------------------------------------------------
#print(log__graph())
hs, vals = log__graph()

#print(hs)
#print(vals)
plt.subplot(212)
plt.xlabel("log(h)")
plt.ylabel("log(err)")
line2 = plt.plot(hs, vals, color="aqua")
plt.grid(True)
plt.show()