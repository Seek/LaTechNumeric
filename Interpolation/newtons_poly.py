import numpy as np
import matplotlib.pyplot as plt
plt.style.use('bmh')

def newtonpoly(x,y, ea):
	assert len(x) == len(y), "x and y must have the same dimensions"
	dim = len(x)
	ddt = np.zeros((dim,dim))
	ddt[:,0] = y
	#Make the ddt
	for i in range(1,dim):
		print("i=" + str(i))
		for j in range (0, dim-i):
			ddt[j,i] = (ddt[j+1,i-1] - ddt[j,i-1])/(x[j+i]-x[j])
	print(ddt)
	
	
x = np.linspace(-2,2, 4)
y = x**3+x**2
print(x)
print(y)
newtonpoly(x, y, 0.1)
#fig, ax = plt.subplots()
#ax.plot(x, y, 'o')
#plt.show()