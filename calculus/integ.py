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
	
print(trap(0.8, 0.2, 0.232))
#Data was taken as f(0), f(0.4), f(0.8); a = 0, b = 0.8, n = 2 or (3-1)
y = [0.2, 2.456, 0.232]
print(trapm(0.4, y))