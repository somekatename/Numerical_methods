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
print()
#Вручную задаем равноотстоящие узлы из интервала непрерывности
nodes = [0]*quantity
values = [0]*quantity
for j in range(quantity):
    nodes[j] = a + j * step
for j in range(quantity):
    values[j] = nodes[j] * m.log(nodes[j] + 1)
print('***Первый способ***')
print(f"Узлы полинома: {nodes}")
print(f"Значения функции: {values}")

#Строим полином Ньютона по заданным вручную узлам:
newton = values[0]
a = 0
mult = (x - nodes[0])
for i in range(1, quantity):
    a = float(((values[i] - newton) / mult).subs({x: nodes[i]}))
    newton += a * mult
    mult *= (x - nodes[i])
print(f"Полином Ньютона для узлов, заданных вручную: N(x) = ")
print(newton)
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
#Строим полином Ньютона по узлам, заданным по формуле:
newton1 = values1[0]
a = 0
mult = (x - nodes1[0])
for i in range(1, quantity):
    a = float(((values1[i] - newton1) / mult).subs({x: nodes1[i]}))
    newton1 += a * mult
    mult *= (x - nodes1[i])
print(f"Полином Ньютона для узлов, заданных по формуле: N(x) = ")
print(newton1)
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

#Строим полином Ньютона (заданный вручную)
plt.subplot(2, 3, 2)
plt.grid()
y2g = [newton.subs(x, i) for i in xg]
plt.title ('Полином Ньютона (первый способ)')
plt.plot(xg, y2g, 'b-')

#Строим полином Ньютона (заданный по формуле)
plt.subplot(2, 3, 3)
plt.grid()
y3g = [newton1.subs(x, i) for i in xg]
plt.title ('Полином Ньютона (второй способ)')
plt.plot(xg, y3g, 'g-')

#Строим полиномы в одном окне
plt.subplot(2, 3, 4)
plt.grid()
plt.plot(xg, y1g, 'r-')
plt.plot(xg, y2g, 'b-')
plt.plot(xg, y3g, 'g-')
plt.legend (('График функции', 'Первый полином Ньютона', 'Второй полином Ньютона'))

#Строим график абсолютной погрешности полинома Ньютона (первый способ)
plt.subplot(2, 3, 5)
plt.grid()
y4g = []
xg1 = []
for i in xg:
    if i * m.log(i + 1) != newton.subs(x, i):
        y4g.append(-m.log10(abs(i * m.log(i + 1) - newton.subs(x, i))))
        xg1.append(i)
plt.title ('График погрешности для первого способа')
plt.plot(xg1, y4g, 'm-')

#Строим график абсолютной погрешности полинома Ньютона (второй способ)
plt.subplot(2, 3, 6)
plt.grid()
y5g = []
xg2 = []
for i in xg:
    if i * m.log(i + 1) != newton1.subs(x, i):
        y5g.append(-m.log10(abs(i * m.log(i + 1) - newton1.subs(x, i))))
        xg2.append(i)
plt.title ('График погрешности для второго способа')
plt.plot(xg2, y5g, 'c-')

plt.show()