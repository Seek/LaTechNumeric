import numpy as np
import pdb
import copy
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
	
def gauss(A, b, tol, err):
	assert A.shape[0] == len(b), "A and b must have the same length"
	dim = A.shape[0]
	x = np.zeros(dim)
	pv = np.arange(0, dim, 1)
	err = 0
	#Eliminate everything but the last row (dim-1)
	for i in range(dim-1):
		#Store the current pivot from the pivot list
		pvt = pv[i]
		#Store the value of the current pivot
		pvv = A[pvt, i]
		#Search the other row specified in the pivot list
		for k in pv:
			#Check if the other rows have larger pivot values
			val = A[k, i]
			#print("val ({0}) > pvv({1})".format(val, pvv))
			if val > pvv:
				#We found a larger row, store the value  and so we can check the others
				pvv = val
				pvt = k
		#Did we find a new pivot?
		if pvt != pv[i]:
			#If we did switch the indices in the pivot list
			#print("We switched row {0} with pivot {1} for row {2} with pivot {3}".format(pv[i], A[pv[i], i], pvt, A[pvt,i]))
			tmp = pv[i]
			pv[i] = pvt
			pv[pvt] = tmp
			#print(pv)
		#Check if the current pivot is close to 0
		#if it is, break and set the error flag
		if 	np.abs(A[pv[i], i]) < tol:
			err = -1
			break
		#Here we actually perform the actual elimination
		for j in range(i+1,dim):
			#print("c = {0}/{1}".format(A[pv[j], i], A[pv[i], i]))
			c = A[pv[j], i]/A[pv[i], i]
			#print(A[pv[j], i:])
			#print((c * A[pv[i], i:]))
			A[pv[j], i:] -= (c * A[pv[i], i:])
			#print(A[pv[j], :])
			b[pv[j]] -= (c*b[pv[i]])
	print(A)
	#Quit here is the system is singular
	if err == -1:
		return x
	#Now we begin back substitution by calculating the last x value
	x[pv[-1]] = b[pv[-1]]/A[pv[-1], -1]
	print(x[pv[-1]])
	#Now we solve the remaining equations
	#dim-2 starts means we begin at second row from the end and go until the 0th row
	for i in range(dim - 2, -1, -1):
		#Grab the corresponding b value
		sum = b[pv[i]]
		#Now we sum from the last column (dim -1 ) to the current column (i-1)
		for j in range(dim - 1, i - 1, -1):
			sum -= x[j] * A[pv[i], j]
		x[i] = sum / A[pv[i], i]
	return x	
e=0
print(gauss1(np.copy(A), np.copy(b)))
print(gauss(np.copy(A), np.copy(b), 0.001, e))
print(b)

	

	