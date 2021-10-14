import numpy as np
import matplotlib.pyplot as plt
import math
from random import uniform
from sympy import *

plt.style.use('classic')
plt.style.use('dark_background')


# f(x) = (x/10)^(sin(x)), x from [0; 10]

def init_func(x):
    return (x / 10)**(math.sin(x))

def inter_func(xs, ys, x):
    res = 0
    for i in range(len(xs)):
        pls = ys[i]
        for j in range(len(xs)):
            if xs[j] != xs[i]:
                c = x - xs[j]
                d = xs[i] - xs[j]
                pls = pls * c / d
        res += pls
    return res    

# Наш исходный график
correct_x = np.linspace(0, 10, 1000)
correct_y = [init_func(i) for i in correct_x]

# Границы отрезка
a = 0
b = 10

h = float(input("Введите ширину h сетки: "))
n = int(10 // h) + 1



# Задаём наше разбиение с использованием исходной функции
first_x = np.linspace(0, 10, num = n)
#print(first_x)
init_x = []
for i in range(n - 1):
    init_x.append((first_x[i] + first_x[i + 1]) / 2)
print("Равномерные узлы: ", end="")
print(init_x)
init_y = []
for i in range(n - 1):
    init_y.append(init_func(init_x[0] + h * i))

#print(init_x)
#print(init_y)

# Получаем точки новой функции через интерполяционный многочлен
out_x = np.linspace(0, 10, 1000)    
out_y1 = [inter_func(init_x, init_y, i) for i in out_x]

#plt.plot(out_x, out_y1, linewidth = 1, color="aqua")
plt.grid(True)

plt.plot(correct_x, correct_y, linewidth = 1, color="yellow")


#----------------------------------------------------------------------------------
# Чебышева узлы на [0, 10]

n = int(input("Введите степень многочлена Чебышева: "))

# Нахождение узлов Чебышева и значений функции в них
cheb_x = []
cheb_y = []
for i in range(1, n + 1):
    cheb_x.append(5 + 5*math.cos((2 * i - 1) * math.pi / (2 * n)))
    cheb_y.append(init_func(cheb_x[i - 1])) 

cheb_x = cheb_x[::-1]
cheb_y = cheb_y[::-1]
print("Узлы Чебышева: ", end="") 
print(cheb_x)
#print(cheb_y)

out_x = np.linspace(0, 10, 1000)
out_y2 = [inter_func(cheb_x, cheb_y, i) for i in out_x]

type = input("Введите желаемый тип: (EQ, CH)  ")

plt.subplot(211)
title = "Равномерные узлы" if type == "EQ" else "Узлы Чебышева"
plt.title(title)

if type == "CH":
    plt.plot(cheb_x, cheb_y, "o", out_x, out_y2, linewidth = 1, color="aqua")
elif type == "EQ":
    plt.plot(init_x, init_y, "o", out_x, out_y1, linewidth = 1, color="aqua")

#plt.plot(out_x, out_y2, linewidth = 1, color="aqua")
plt.grid(True)

plt.plot(correct_x, correct_y, linewidth = 1, color="yellow")

plt.subplot(212)
plt.title("Ошибка")
if type == "EQ":
    plt.plot(correct_x, [correct_y[i] - out_y1[i] for i in range(len(out_y1))], linewidth = 1, color="aqua")
elif type == "CH":
    plt.plot(correct_x, [correct_y[i] - out_y2[i] for i in range(len(out_y2))], linewidth = 1, color="yellow")

plt.grid(True)

plt.show()    



#---------------------------------------------------------------------------------
# Функция ошибки
"""
n = int(input("Введите желаемый порядок дифференцирования: "))
x = Symbol('x')
func = (x / 10)**(sin(x))

#er_x = np.linspace(0, 10, n)
#er_y = [init_func(i) for i in er_x]

def added_diff(xs, n, f = func):
    res = 1
    for i in range(n):
        f = f.diff()
    f1 = lambdify(x, f, 'numpy')    
    ksi = uniform(xs[n - 1], xs[n])
    res = res * f1(ksi) / math.factorial(n)
    return res        
        
def diff_res(xs, er_xs, n):
    res = 0
    for i in range(n):
        res += added_diff(xs, i)
    return res

#for i in range(n):
#    func = func.diff(x)

#print(func)
#ksi = uniform(0, 10)
#print(ksi)

#f = lambdify(x, func, 'numpy')



er_x = np.linspace(0, 10, 100)


print("", end="")  
"""

