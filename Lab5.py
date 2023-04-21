import math as m
import numpy as np
import matplotlib.pyplot as plt
from sympy import *

# Задаем нашу функцию y = xln(x+1):
x = Symbol('x', real=True)
y = x * log(x + 1)
print("Введите начало интервала непрерывности")
a = float(input())
print("Введите конец интервала непрерывности")
b = float(input())
print("Введите расстояние между равноотстоящими узлами")
step = float(input())

# количество узлов зависит от расстояния между ними
quantity = round((b - a) / step + 1)
print(f"Количество узлов интерполирования: {quantity}")

#Вручную задаем равноотстоящие узлы из интервала непрерывности [0, 4.5]
nodes = [0]*quantity
values = [0]*quantity
for j in range(quantity):
    nodes[j] = a + j * step
for j in range(quantity):
    values[j] = nodes[j] * m.log(nodes[j] + 1)
print('***Первый способ***')
print(f"Узлы полинома: {nodes}")
print(f"Значения функции: {values}")


#Строим полином Лагранжа по заданным вручную узлам:
lagrange = 0
index = 0
numerator = 1
denominator = 1
for i in range(quantity):
    for j in range(quantity):
        if j == index:
            continue
        numerator *= (x - nodes[j])
        denominator *= (nodes[index] - nodes[j])
    numerator *= values[i]
    lagrange += numerator / denominator
    numerator = 1
    denominator = 1
    index += 1
print(f"Полином Лагранжа для узлов, заданных вручную: L(x) = {lagrange}")
print()

#Задаем таблицу для построения полинома по формуле 3.2:
nodes1 = [0]*quantity
values1 = [0]*quantity
for j in range(quantity):
    nodes1[j] = 0.5 * ((b - a) * m.cos((2 * j + 1) * m.pi
                        / (2 * (quantity + 1))) + (b + a))
for j in range(quantity):
    values1[j] = nodes1[j] * m.log(nodes1[j] + 1)
print('***Второй способ***')
print(f"Узлы полинома: {nodes1}")
print(f"Значения функции: {values1}")

#Строим полином Лагранжа по узлам, заданным по формуле:
lagrange1 = 0
index = 0
numerator = 1
denominator = 1
for i in range(quantity):
    for j in range(quantity):
        if j == index:
            continue
        numerator *= (x - nodes1[j])
        denominator *= (nodes1[index] - nodes1[j])
    numerator *= values1[i]
    lagrange1 += numerator / denominator
    numerator = 1
    denominator = 1
    index += 1
print(f"Полином Лагранжа для узлов, заданных по формуле: L(x) = {lagrange1}")
print()

# ***ГРАФИЧЕСКОЕ ПРЕДСТАВЛЕНИЕ РЕЗУЛЬТАТОВ***
xg = np.linspace(a, b, 100)
plt.axis([a, b, 0, b*m.log(b + 1)])

#Строим график нашей функции
plt.subplot(2, 3, 1)
plt.grid()
y1g = [i * m.log(i + 1) for i in xg]
plt.title ('График функции')
plt.plot(xg, y1g, 'r-')

#Строим полином Лагранжа (заданный вручную)
plt.subplot(2, 3, 2)
plt.grid()
y2g = [lagrange.subs(x, i) for i in xg]
plt.title ('Полином Лагранжа (первый способ)')
plt.plot(xg, y2g, 'b-')

#Строим полином Лагранжа (заданный по формуле)
plt.subplot(2, 3, 3)
plt.grid()
y3g = [lagrange1.subs(x, i) for i in xg]
plt.title ('Полином Лагранжа (второй способ)')
plt.plot(xg, y3g, 'g-')

#Строим полиномы в одном окне
plt.subplot(2, 3, 4)
plt.grid()
plt.plot(xg, y1g, 'r-')
plt.plot(xg, y2g, 'b-')
plt.plot(xg, y3g, 'g-')
plt.legend (('График функции', 'Первый полином Лагранжа', 'Второй полином Лагранжа'))

#Строим график абсолютной погрешности полинома Лагранжа (первый способ)
plt.subplot(2, 3, 5)
plt.grid()
y4g = []
xg1 = []
for i in xg:
    if i * m.log(i + 1) != lagrange.subs(x, i):
        y4g.append(-m.log10(abs(i * m.log(i + 1) - lagrange.subs(x, i))))
        xg1.append(i)
plt.title ('График погрешности для первого способа')
plt.plot(xg1, y4g, 'm-')

#Строим график абсолютной погрешности полинома Лагранжа (второй способ)
plt.subplot(2, 3, 6)
plt.grid()
y5g = []
xg2 = []
for i in xg:
    if i * m.log(i + 1) != lagrange1.subs(x, i):
        y5g.append(-m.log10(abs(i * m.log(i + 1) - lagrange1.subs(x, i))))
        xg2.append(i)
plt.title ('График погрешности для второго способа')
plt.plot(xg2, y5g, 'c-')

plt.show()