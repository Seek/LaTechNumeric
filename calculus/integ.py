import numpy as np
def trap(h, y1, y2):
	""" h is defined as (b-a) in which a and b are the limits of integration """
	return h * (y1+y2)/2 

def trapm(h, y):
	""" h = (b-a)/(len(y)-1) where b and a are the limits of integration
		y  is an numpy compatability array
	"""
	n = len(y)-1
	sum = y[0] + y[-1]
	sum += np.sum(y[1:-1]) * 2
	return h * (sum/2)
	
def simp13(h, y1, y2, y3):
	""" h is defined as (b-a) in which a and b are the limits of integration """
	return (h) * (y1 + 4*y2 + y3)/6

def simp13m(h, y):
	""" h is defined as (b-a)/(len(y)-1) in which a and b are the limits of integration """
	sum = y[0] + y[-1]
	sum += np.sum(y[1:-1:2]) * 4
	sum += np.sum(y[2:-1:2]) * 2
	return sum * h / 3
	
def simp38(h, y1, y2, y3, y4):
	""" h is defined as (b-a) in which a and b are the limits of integration """
	return h * (y1 + 3*y2 + 3*y3 + y4)/8

	
print(trap(0.8, 0.2, 0.232))
#Data was taken as f(0), f(0.4), f(0.8); a = 0, b = 0.8, n = 2 or (3-1)
y = [0.2, 2.456, 0.232]
print(trapm(0.4, y))
print(simp13(0.8, *y))
#Data was taken as f(0), f(0.4), f(0.8); a = 0, b = 0.8, n = 2 or (3-1)
y2 = [0.2, 1.288, 2.456, 3.464, 0.232]
print(simp13m(0.8/4, y2))
y3 = [0.2, 1.432724, 3.487177, 0.232]
print(simp38(0.8, *y3))