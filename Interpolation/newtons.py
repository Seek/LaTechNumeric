import numpy as np

def newtonpoly(x,y, xrange):
	assert len(x) == len(y), "x and y must have the same dimensions"
	dim = len(x)
	ddt = np.zeros((dim,dim))
	ddt[:,0] = y
	#Make the ddt
	for i in range(1,dim):
		for j in range (0, dim-i):
			ddt[j,i] = (ddt[j+1,i-1] - ddt[j,i-1])/(x[j+i]-x[j])
	#Evaluate the polynomial
	yint = np.zeros(len(xrange))
	for i in range(dim):
		tx = np.ones(len(xrange))
		for j in range(i):
			tx *= (xrange-x[j])
		tx *= ddt[0,i]
		yint += tx
	return yint