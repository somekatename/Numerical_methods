# g(x) = ctg(x)
# h(x) = x^2
import math as m
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from sympy import *

print('****Локализация корня методом Ньютона****')

x = Symbol('x', real=True)
# Уравнение и его производная
f = tan(0.5*x+0.2) - x**2
df = 1/(2*(cos(0.5*x+0.2))**2) - 2*x
#Начальное приближение, полученное графическим методом
x0 = -0.1
x1 = x0 - float(f.subs({x: x0}))/ float(df.subs({x: x0}))
while abs(x1 - x0) > 0.0001:
    x0 = x1
    x1 = x0 - float(f.subs({x: x0}))/ float(df.subs({x: x0}))
print(x1)
#Построение графика
xg = np.linspace(-5, 5, 100)
y1g = [i ** 2 for i in xg]
y2g = [m.tan(0.5*j+0.2) for j in xg]

plt.axis([-5, 5, 0, 10])
plt.xlabel('x') # ось абсцисс
plt.ylabel('y1, y2') # ось ординат
plt.grid() # включение отображение сетки
plt.plot(xg, y1g, 'r-', xg, y2g, 'b-')
plt.show()



