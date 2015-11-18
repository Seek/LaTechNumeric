import numpy as np
import pdb
A = np.array([
[1,1,-2,1,3,-1],
[2,-1,1,2,1,-3],
[1,3,-3,-1,2,1],
[5,2,-1,-1,2,1],
[-3,-1,2,3,1,3],
[4,3,1,-6,-3,-2]
], dtype=float)
AB = np.array([
[1,1,-2,1,3,-1],
[2,-1,1,2,1,-3],
[1,3,-3,-1,2,1],
[5,2,-1,-1,2,1],
[-3,-1,2,3,1,3],
[4,3,1,-6,-3,-2],
[4,20,-15,-3, 16, -27]
], dtype=float)
b = np.array([4,20,-15,-3, 16, -27], dtype=float)

#Simplest version
def gauss1(A, b):
	assert A.shape[0] == len(b), "A and b must have the same length"
	dim = A.shape[0]
	x = np.zeros(dim)
	#Elimination
	for i in range(dim-1):
		for j in range(i+1,dim):
			c = A[j, i]/A[i, i]
			A[j, :] -= (c*A[i, :])
			b[j] -= (c*b[i])
	#Substitution
	x[-1] = b[-1]/A[-1, -1]
	for i in range(dim - 2, -1, -1):
		sum = b[i]
		for j in range(dim - 1, i - 1, -1):
			sum -= x[j] * A[i, j]
		x[i] = sum / A[i, i]
	return x

#With pivoting, we make copies though.
def gauss2(A, b, tol, err):
	assert A.shape[0] == len(b), "A and b must have the same length"
	dim = A.shape[0]
	x = np.zeros(dim)
	#Elimination
	for i in range(dim-1):
		#Partial pivoting
		p = np.argmax(A[:, i])
		if p < i:
			print("Switching row {0} with {1}".format(i, p))
			print("Row {0} has pivot {1} the current pivot is {2}".format(p, A[p, i], A[i,i]))
			tmp = np.copy(A[i, :])
			A[i, :] = A[p, :]
			A[p, :] = np.copy(tmp)
			tmp = b[i]
			b[i] = b[p]
			b[p] = tmp
		#Actual elimination
		if 	np.abs(A[i,i]) < tol:
			err = -1
			break
		for j in range(i+1,dim):
			c = A[j, i]/A[i, i]
			A[j, :] -= (c*A[i, :])
			b[j] -= (c*b[i])
	#Substitution
	if err == -1:
		return x
	x[-1] = b[-1]/A[-1, -1]
	for i in range(dim - 2, -1, -1):
		sum = b[i]
		for j in range(dim - 1, i - 1, -1):
			sum -= x[j] * A[i, j]
		x[i] = sum / A[i, i]
	return x

def gaussjordan(A, tol, err):
	dim = A.shape[0]
	#Elimination
	for i in range(dim-1):
		#Partial pivoting
		p = np.argmax(A[:, i])
		if p < i:
			print("Switching row {0} with {1}".format(i, p))
			print("Row {0} has pivot {1} the current pivot is {2}".format(p, A[p, i], A[i,i]))
			tmp = np.copy(A[i, :])
			A[i, :] = A[p, :]
			A[p, :] = np.copy(tmp)
			tmp = b[i]
			b[i] = b[p]
			b[p] = tmp
		#Actual elimination
		if 	np.abs(A[i,i]) < tol:
			err = -1
			break
		for j in range(i+1,dim):
			A[i, :] /= A[i,i]
			c = A[j, i]/A[i, i]
			A[j, :] -= (c*A[i, :])
	A[-1, :] /= A[-1,-1]
	print(A)
	#Substitution
	return A
e = 0
print(AB)
#x =  gauss1(A, b)
x = gauss1(np.copy(A), np.copy(b))
x1 = gauss2(np.copy(A), np.copy(b), 0.001, e)
aa = gaussjordan(AB, 0.001, e)
print(b)
print(np.dot(A, x))
print(np.dot(A, x1))

	

	