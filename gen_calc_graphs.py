import numpy as np
import matplotlib.pyplot as plt
import calculus.integ as numint
import scipy.integrate as spint
from interpolation.lagrange import lagrangepoly
#Real value 7.32473
def f(x):
	return -x**5 + 5*x**4 - 7*x**3 + x**2 + 4*x + 0.5

x = np.linspace(-5,5, 150)
y = f(x)
a = 0.5
b = 2.6
x2 = np.linspace(a, b)
fig, ax = plt.subplots()
ax.plot(x, y)
ax.set_xlim(-3, 5)
ax.set_ylim(-4, 7)
ax.axhline(0, color='black')
ax.legend(['f(x)'])
ax.fill_between(x2, 0, f(x2), color = 'g', alpha=0.3)
plt.title('Integral of f(x)')

#Trap rule
fig, ax = plt.subplots()
ax.plot(x, y)
ax.set_xlim(-3, 5)
ax.set_ylim(-4, 7)
ax.axhline(0, color='black')
ax.legend(['f(x)'])

h = b-a

y2 = f(a) + (f(b) - f(a))/(b - a) * (x2 - a)
ax.plot(x2, y2, '-g')
ax.fill_between(x2, 0, y2, color = 'g', alpha=0.3)
ax.grid(True)
textstr = str(numint.trap(h, f(a), f(b)))
props = dict(boxstyle='round', facecolor='white')
# place a text box in upper left in axes coords
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

plt.title('Trapezoid Rule')

#Trap rule with multiple segs
fig, ax = plt.subplots()
ax.plot(x, y)
numpts = 5
n = 5-1
ax.set_xlim(-3, 5)
ax.set_ylim(-4, 7)
ax.axhline(0, color='black')
ax.legend(['f(x)'])
x3 = np.linspace(a, b, numpts)
x4 = np.linspace(a, b, numpts)
x3 = x3[1:]
lastx = a
for i in x3:
	xn = np.linspace(lastx, i, 15)
	yn = f(lastx) + (f(i) - f(lastx))/(i - lastx) * (xn - lastx)
	ax.plot(xn, yn, '-g')
	ax.fill_between(xn, 0, yn, color='g', alpha=0.3)
	lastx = i
	
textstr = str(numint.trapm(h/n, f(x4)))
props = dict(boxstyle='round', facecolor='white')
# place a text box in upper left in axes coords
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)
ax.grid(True)
plt.title('Trapezoid Rule w/ Multiple Segments')

#Simps 1/3 rule
fig, ax = plt.subplots()
ax.plot(x, y)
ax.set_xlim(-3, 5)
ax.set_ylim(-4, 7)
ax.axhline(0, color='black')
ax.legend(['f(x)'])
a = 0.5
b = 2.6
h = b-a
xx = np.linspace(a, b, 3)
x2 = np.linspace(a, b)
y2 = lagrangepoly(xx, f(xx), x2)
ax.plot(x2, y2, '-g')
ax.fill_between(x2, 0, y2, color = 'g', alpha=0.3)
ax.grid(True)
textstr = str(numint.simp13(h, f(xx[0]), f(xx[1]), f(xx[2])))
props = dict(boxstyle='round', facecolor='white')
# place a text box in upper left in axes coords
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

plt.title('Simpson 1/3 Rule')

#Simps 1/3 rule with multiple segs
fig, ax = plt.subplots()
ax.plot(x, y)
numpts = 5
n = 5-1
ax.set_xlim(-3, 5)
ax.set_ylim(-4, 7)
ax.axhline(0, color='black')
ax.legend(['f(x)'])
x3 = np.linspace(a, b, numpts)
x4 = np.linspace(a, b, numpts)
x3 = x3[1:]
lastx = a
for i in x3:
	xn = np.linspace(lastx, i, 15)
	xx2 = np.linspace(lastx, i, 3)
	yn = lagrangepoly(xx2, f(xx2), xn)
	ax.plot(xn, yn, '-g')
	ax.fill_between(xn, 0, yn, color='g', alpha=0.3)
	lastx = i
	
textstr = str(numint.simp13m(h/n, f(x4)))
props = dict(boxstyle='round', facecolor='white')
# place a text box in upper left in axes coords
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)
ax.grid(True)
plt.title('Simpson Rule 1/3 w/ Multiple Segments')
plt.show()
