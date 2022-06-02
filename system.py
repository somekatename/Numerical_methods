import math as m
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from sympy import *

x = Symbol('x', real=True)
y = Symbol('y', real=True)

x0 = -0.8
y0 = -0.1
eps = 0.0001
# Выписываем уравнения системы и элементы якобиана
f1 = sin(y + 0.5) - x - 1
f2 = y + cos(x-2)
f1x = -1.0
f1y = cos(y+0.5)
f2x = -sin(x-2)
f2y = 1.0

# Находим Якобиан в точке начального приближения и решаем систему
a12 = f1y.subs({y: y0})
a21 = f2x.subs({x: x0})
b1 = -f1.subs({x: x0, y: y0})
b2 = -f2.subs({x: x0, y: y0})
A = np.array([[f1x, a12], [a21, f2y]], dtype=np.float64)
B = np.array([b1, b2], dtype=np.float64)
C = np.linalg.solve(A, B)
x1 = x0 + C[0]
y1 = y0 + C[1]

while sqrt((x1 - x0)**2 + (y1 - y0)**2) > eps:
    x0 = x1
    y0 = y1
    a12 = f1y.subs({y: y0})
    a21 = f2x.subs({x: x0})
    b1 = -f1.subs({x: x0, y: y0})
    b2 = -f2.subs({x: x0, y: y0})
    A = np.array([[f1x, a12], [a21, f2y]], dtype=np.float64)
    B = np.array([b1, b2], dtype=np.float64)
    C = np.linalg.solve(A, B)
    x1 = x0 + C[0]
    y1 = y0 + C[1]
print('x =', x1, '\n', 'y =', y1)

xg = np.linspace(-1.3, 0, 100)
y1g = [-m.cos(i-2) for i in xg]
y2g = [m.asin(1 + j) - 0.5 for j in xg]
plt.axis([-1.25, 0, -1, 1])
plt.xlabel('x')
plt.ylabel('y1, y2')
plt.grid()
plt.plot(xg, y1g, 'r-', xg, y2g, 'b-')
plt.show()


