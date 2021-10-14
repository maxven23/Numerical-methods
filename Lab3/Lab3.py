import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import *

#-------------------------------------------------------

plt.style.use('classic')
plt.style.use('dark_background')

#-------------------------------------------------------

a = -2
b = 1

#-------------------------------------------------------

def init_func(x):
    return x / sqrt(5 - x*x)
#def init_func(x):
#    return x*(math.e**x)

def rect_method(n, func=init_func):
    h = (b - a) / n
    X = np.linspace(a, b, n)
    res = 0
    for i in range(1, n):
        res += init_func(X[i] - h / 2)
    return res * h

def trapez_method(n, func=init_func):
    h = (b - a) / n
    X = np.linspace(a, b, n)
    Y = [init_func(i) for i in X]
    res = (Y[0] + Y[n - 1]) / 2
    for i in range(1, n - 1):
        res += init_func(X[i])
    return res * h

def Simpson_method(n, func=init_func):
    h = (b - a) / (2 * n)
    X = np.linspace(a, b, (2 * n + 1))
    Y = [init_func(i) for i in X]
    res = Y[0] + Y[2 * n]
    for i in range(1, 2 * n + 1, 2):
        res += 4*Y[i]
    for j in range(2, 2 * n, 2):
        res += 2*Y[j]
    return (res * h / 3)

#-------------------------------------------------------

correct = -1
#correct = 3 / (math.e * math.e)

#-------------------------------------------------------

n = int(input("Введите количество делений сетки: "))
print()
print("Точное значение: \t", correct)
print("Метод прямоугольников: \t", rect_method(n))
print("Метод трапеций: \t", trapez_method(n))
print("Метод Симпсона: \t", Simpson_method(n))

#-------------------------------------------------------

def rect_graph(n, func=init_func):
    h = (b - a) / (n - 1)
    X = np.linspace(a, b, n)

    Y = []
    Y.append(init_func(a))
    
    F = [init_func(X[i] - h/2) for i in range(1, n)]
    Y.extend(F)

    MINS = [0 for i in X]

    plt.step(X, Y, where='pre', color="white")
    plt.title("Метод прямоугольников")
    plt.grid(True)
    
    plt.vlines(X, MINS, Y, color="white", linewidth=1)
    plt.axhline(y=0, xmin=X[0], xmax=X[len(X) - 1], color="white", linewidth=2)
    
    Xs = np.linspace(a, b, num=1000)
    Ys = [init_func(i) for i in Xs]
    
    plt.plot(Xs, Ys, color="lime", linewidth=1, linestyle="--", label="Исходная функция")
    plt.legend(loc='lower right')
    plt.fill_between(X, np.array(Y, dtype="float"), 0, step="pre", color="white", alpha=0.25)
    plt.show()

#-------------------------------------------------------
    
def parabola(a, b, c, x):
    return a*x*x + b * x + c
def Simpson_graph(n, func=init_func):
    h = (b - a) / n
    
    X = np.linspace(a, b, n)

    Y = [init_func(i) for i in X]
#    print(X)
#    print(Y)
    MINS = [0 for i in X]
    
    for i in range(1, n):
        tempX = []
        tempX.append(X[i - 1])
        tempX.append((X[i - 1] + X[i])/2)
        tempX.append(X[i])
        tempY = [init_func(i) for i in tempX]
        A = (tempY[2] - (tempX[2]*(tempY[1] - tempY[0]) + tempX[1]*tempY[0] - tempX[0]*tempY[1]) / (tempX[1] - tempX[0])) / (tempX[2] * (tempX[2] - tempX[1] - tempX[0]) + tempX[1]*tempX[0])
        B = (tempY[1] - tempY[0]) / (tempX[1] - tempX[0]) - A * (tempX[1] + tempX[0])
        C = (tempX[1]*tempY[0] - tempX[0]*tempY[1]) / (tempX[1] - tempX[0]) + A * tempX[0] * tempX[1]
        goodX = np.linspace(X[i - 1], X[i], num=100)
        goodY = [parabola(A, B, C, i) for i in goodX]
        plt.plot(goodX, goodY, color="white", linewidth=1)
        plt.fill_between(goodX, np.array(goodY, dtype="float"), color="white", alpha=0.25)
        

    plt.title("Метод Симпсона")
    plt.grid(True)
    
    plt.vlines(X, MINS, Y, color="white", linewidth=1)
    plt.axhline(y=0, xmin=X[0], xmax=X[len(X) - 1], color="white", linewidth=2)
    
    Xs = np.linspace(a, b, num=1000)
    Ys = [init_func(i) for i in Xs]
    
    plt.plot(Xs, Ys, color="lime", linewidth=1, linestyle="--", label="Исходная функция")
    plt.legend(loc='lower right')
    plt.show()

#-------------------------------------------------------
    
def trapez_graph(n, func=init_func):
    h = (b - a) / n
    
    X = np.linspace(a, b, n)
    
    Y = []
    
    Y = [init_func(i) for i in X]
#    print(X)
#    print(Y)
    MINS = [0 for i in X]

    plt.plot(X, Y, color="white")
    plt.title("Метод трапеций")
    plt.grid(True)
    
    plt.vlines(X, MINS, Y, color="white", linewidth=1)
    plt.axhline(y=0, xmin=X[0], xmax=X[len(X) - 1], color="white", linewidth=2)
    
    Xs = np.linspace(a, b, num=1000)
    Ys = [init_func(i) for i in Xs]
    
    plt.plot(Xs, Ys, color="lime", linewidth=1, linestyle="--", label="Исходная функция")
    plt.legend(loc='lower right')   
    plt.fill_between(X, np.array(Y, dtype="float"), 0, color="white", alpha=0.25)
    plt.show()
 
#-------------------------------------------------------   

rect_graph(n)
input("Enter any key: ")
trapez_graph(n)
input("Enter any key: ")
Simpson_graph(n)
input("Enter any key: ")

hs = []
rectY = []
trapezY = []
simpsonY = []
for i in range(3, 1000, 2):
    hs.append(((b - a) / i))
    rectY.append(log(abs(rect_method(i) - correct)))
    trapezY.append(log(abs(trapez_method(i) - correct)))
    simpsonY.append(log(abs(Simpson_method(i) - correct)))
       
plt.title("Зависимость ошибки от размера шага")
plt.xlabel("h")
plt.ylabel("log(err)")
plt.grid(True)
plt.plot(hs, rectY, color="white", linewidth=2, label="Rectangles method")
plt.plot(hs, trapezY, color="aqua", linewidth=2, label="Trapezes method")
plt.plot(hs, simpsonY, color="lime", linewidth=2, label="Simpson method")
plt.legend(loc='lower right')
plt.show()