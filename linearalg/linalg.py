import numpy as np

A = np.array([[3, -0.1, -0.2], [0.1, 7, -0.3], [0.3, -0.2, 10]])
b = [7.85, -19.3, 71.4]
def gauss1(A, b):
	assert A.shape[0] == len(b), "A and b must have the same length"
	dim = A.shape[0]
	x = np.zeros(dim)
	tx = 0
	for i in range(dim-1):
		for j in range(i+1,dim):
			c = A[j, i]/A[i, i]
			A[j, :] -= (c*A[i, :])
			b[j] -= (c*b[i])
	x[-1] = b[-1]/A[-1, -1]
	for i in range(dim - 2, -1, -1):
		print("i = " + str(i))
		sum = b[i]
		for j in range(dim - 1, i - 1, -1):
			sum -= x[j] * A[i, j]
		x[i] = sum / A[i, i]
	return x

x = gauss1(A, b)
print(x)
	

	