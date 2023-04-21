import math as m
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from sympy import *

def transpose_matrix(lst):
    matrix = list(map(list, zip(*lst)))
    return matrix

def matrix_multiplication(lst1, lst2):
    n1, m1 = len(lst1), len(lst1[0])
    n2, m2 = len(lst2), len(lst2[0])
    matrix_mult = [[0] * m2 for _ in range(n1)]
    for i in range(n1):
        for j in range(m2):
            for s in range(n2):
                matrix_mult[i][j] += float(lst1[i][s] * lst2[s][j])
    return matrix_mult

def printing_matrix(matrix):
    for row in matrix:
        for elem in row:
            print(f"{elem:.{3}f}", end=' ')
        print()

# Задаем нашу функцию y = xln(x+1):
x = Symbol('x', real=True)
a0 = Symbol('a0', real=True)
a1 = Symbol('a1', real=True)
a2 = Symbol('a2', real=True)
a3 = Symbol('a3', real=True)
y = x * log(x + 2)
print("Введите начало интервала непрерывности")
start = float(input())
print("Введите конец интервала непрерывности")
end = float(input())
print("Введите расстояние между равноотстоящими узлами")
step = float(input())

# количество узлов зависит от расстояния между ними
n = round((end - start) / step + 1)
print(f"Количество узлов аппроксимации: {n}")
#------------------------------------------------------------------------------------------------------------#
#                                 МЕТОД НАИМЕНЬШИХ КВАДРАТОВ                                                 #
#------------------------------------------------------------------------------------------------------------#
#Вручную задаем равноотстоящие узлы из интервала непрерывности [-1, 1]
nodes = [0]*n
values = [[0] for _ in range(n)]
for j in range(n):
    nodes[j] = start + j * step

# Для значений функции введем дополнительно погрешность и введем у как вектор-столбец
error = 0.00001
for j in range(n):
    values[j][0] = nodes[j] * m.log(nodes[j] + 2) + error
    error += 0.000001
print(f"Узлы полинома: {nodes}")
print("Значения функции:")
printing_matrix(values)

#Создаем вектор-столбец phi, обобщенный полином Pm:
phi = [[x**0], [x], [x**2], [x**3]]
a = [[a0], [a1], [a2], [a3]]
Pm = 0
for i in range(4):
    Pm += a[i][0] * phi[i][0]
print("Обобщенный полином Pm = ", Pm)

#Введем среднеквадратичное отклонение Pm от y:
sgm = 0
for i in range(n):
    sgm += (Pm.subs({x: nodes[i]}) - values[i][0])**2
q = [[0]*4 for _ in range(n)]
for i in range(n):
    for j in range(4):
        q[i][j] = phi[j][0].subs({x: nodes[i]})
print("Матрица Q:")
printing_matrix(q)
print("Транспонированная Q")
qt = transpose_matrix(q)
printing_matrix(qt)

#Введем матрицу H и вектор b, записываем решение системы Ha = b в переменную a_solved:
h = np.array(matrix_multiplication(qt, q))
b = np.array(sum(matrix_multiplication(qt, values), []))
a_solved = np.linalg.solve(h, b)

#Выводим найденный полином Pm:
Pm = Pm.subs({a0: a_solved[0], a1:a_solved[1], a2:a_solved[2], a3:a_solved[3]})
print(f"Аппроксимирующий полином Pm: {Pm}")
print(f"Величина среднеквадратичного отклонения: {m.sqrt(sgm.subs({a0: a_solved[0], a1: a_solved[1], a2: a_solved[2], a3: a_solved[3]}) / (n + 1))}")

#------------------------------------------------------------------------------------------------------------#
#                                 ПОЛИНОМЫ ЛЕЖАНДРА                                                          #
#------------------------------------------------------------------------------------------------------------#
print()
print("ПОЛИНОМЫ ЛЕЖАНДРА")
#Введем многочлен наилучшего приближения Q3:
Q3 = 0
delta = sym.integrate(y**2, (x, -1, 1))
for i in range(4):
    divider = (m.factorial(i) * 2**i)
    numerator = sym.diff((1 - x**2)**i, x, i)
    c = (sym.integrate(y * numerator / divider, (x, -1, 1))) \
        / (sym.integrate((numerator / divider)**2, (x, -1, 1)))
    Q3 += c * numerator / divider
    delta -= c**2 * (sym.integrate((numerator / divider)**2, (x, -1, 1)))
print("Полином наилучшего приближения Q3 = ", Q3)
print(f"Величина среднеквадратичного отклонения: {m.sqrt(delta)}")


# ***ГРАФИЧЕСКОЕ ПРЕДСТАВЛЕНИЕ РЕЗУЛЬТАТОВ***
xg = np.linspace(start, end, 100)
plt.axis([start, end, 0, end*m.log(end + 2)])

#Строим график нашей функции
plt.subplot(2, 3, 1)
plt.grid()
y1g = [i * m.log(i + 2) for i in xg]
plt.title ('График функции')
plt.plot(xg, y1g, 'r-')

#График аппроксимирующего полинома
plt.subplot(2, 3, 2)
plt.grid()
y2g = [Pm.subs({x: i}) for i in xg]
plt.title ('График аппроксимирующего полинома')
plt.plot(xg, y2g, 'b-')

#График полинома наилучшего приближения с полиномами Лежандра
plt.subplot(2, 3, 3)
plt.grid()
y3g = [Q3.subs({x: i}) for i in xg]
plt.title ('График полинома наилучшего приближения')
plt.plot(xg, y3g, 'g-')

#Графики вместе
plt.subplot(2, 3, 4)
plt.grid()
plt.plot(xg, y1g, 'r-')
plt.plot(xg, y2g, 'b-')
plt.plot(xg, y3g, 'g-')
plt.legend (('График функции', 'Аппроксимирующий полином 1', "Полином наилучшего приближения (Лежандр)"))

#Погрешности каждого из способов:
plt.subplot(2, 3, 5)
plt.grid()
y4g = [-m.log10(abs(i * m.log(i + 2) - Pm.subs({x: i}))) for i in xg]
plt.title ('График погрешности для первого способа')
plt.plot(xg, y4g, 'c-')

plt.subplot(2, 3, 6)
plt.grid()
y5g = [-m.log10(abs(i * m.log(i + 2) - Q3.subs({x: i}))) for i in xg]
plt.title ('График погрешности для второго способа')
plt.plot(xg, y5g, 'm-')

plt.show()