import numpy as np

def lagrangepoly(x, y, xrange):
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