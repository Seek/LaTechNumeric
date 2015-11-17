import numpy as np
import matplotlib.pyplot as plt
plt.style.use('bmh')

def lagrangepoly(x,y, xrange):
	assert len(x) == len(y), "x and y must have the same dimensions"
	dim = len(x)
	dim2 = len(xrange)
	yint = np.zeros(dim2)
	for i in range(dim):
		tt = np.ones(dim2)
		for j in range(dim):
			if i == j:
				continue
			tt *= (xrange - x[j])/(x[i]-x[j])
		yint += y[i] * tt
	return yint
	
x = np.linspace(-np.pi, np.pi, 5)
x2 = np.linspace(-np.pi, np.pi)
y = np.sin(x)
yint = lagrangepoly(x, y, x2)
fig, ax = plt.subplots()
ax.plot(x, y, 'o')
ax.plot(x2, yint, '-')
plt.show()