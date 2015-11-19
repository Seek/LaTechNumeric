import numpy as np
import pdb
from scipy import linalg as splinalg
# A = np.array([
    # [1, 1, -2, 1, 3, -1],
    # [2, -1, 1, 2, 1, -3],
    # [1, 3, -3, -1, 2, 1],
    # [5, 2, -1, -1, 2, 1],
    # [-3, -1, 2, 3, 1, 3],
    # [4, 3, 1, -6, -3, -2]
# ], dtype=float)
# b = np.array([4, 20, -15, -3, 16, -27], dtype=float)

A = np.array([
[8,4,4],
[2,-4,1],
[2,-1,3]
], dtype = float)
b = np.array([
80, 7, 22
], dtype=float)

# A = np.array([
# [3,-0.1,-0.2],
# [0.1,7,-0.3],
# [0.3,-0.2,10]
# ], dtype = float)
# b = np.array([
# 7.85, -19.3, 71.4
# ], dtype=float)



# Simplest version
def gauss1(A, b):
    assert A.shape[0] == len(b), "A and b must have the same length"
    dim = A.shape[0]
    x = np.zeros(dim)
    # Elimination
    for i in range(dim - 1):
        for j in range(i + 1, dim):
            c = A[j, i] / A[i, i]
            A[j, :] -= (c * A[i, :])
            b[j] -= (c * b[i])
    # Substitution
    x[-1] = b[-1] / A[-1, -1]
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
    # Eliminate everything but the last row (dim-1)
    for i in range(dim - 1):
        # Store the current pivot from the pivot list
        pvt = pv[i]
        # Store the value of the current pivot
        pvv = A[pvt, i]
        # Search the other row specified in the pivot list
        for k in pv:
            # Check if the other rows have larger pivot values
            val = A[k, i]
            # print("val ({0}) > pvv({1})".format(val, pvv))
            if val > pvv:
                # We found a larger row, store the value  and so we can check the others
                pvv = val
                pvt = k
        # Did we find a new pivot that is in a row below us?
        if pvt > pv[i]:
            # If we did switch the indices in the pivot list
            #print("We switched row {0} with pivot {1} for row {2} with pivot {3}".format(pv[i], A[pv[i], i], pvt, A[pvt,i]))
            tmp = pv[i]
            pv[i] = pvt
            pv[pvt] = tmp
        # print(pv)
        # Check if the current pivot is close to 0
        # if it is, break and set the error flag
        if np.abs(A[pv[i], i]) < tol:
            err = -1
            break
        # Here we actually perform the actual elimination
        for j in range(i + 1, dim):
            # print("c = {0}/{1}".format(A[pv[j], i], A[pv[i], i]))
            c = A[pv[j], i] / A[pv[i], i]
            # print(A[pv[j], i:])
            # print((c * A[pv[i], i:]))
            A[pv[j], i:] -= (c * A[pv[i], i:])
            # print(A[pv[j], :])
            b[pv[j]] -= (c * b[pv[i]])
    # print(A)
    #print(b)
    # Quit here is the system is singular
    if err == -1:
        return x
    # Now we begin back substitution by calculating the last x value
    x[-1] = b[pv[-1]] / A[pv[-1], -1]
    # Now we solve the remaining equations
    # dim-2 starts means we begin at second row from the end and go until the 0th row
    for i in range(dim - 2, -1, -1):
        # Grab the corresponding b value
        sum = b[pv[i]]
        # Now we sum from the last column (dim -1 ) to the current column (i-1)
        for j in range(dim - 1, i - 1, -1):
            sum -= x[j] * A[pv[i], j]
        x[i] = sum / A[pv[i], i]
    return x

def lu_factor(A, tol, err):
    """Returns the matrix A with the LU matrices and a pivot vector containing information on how the matrix was eliminated.
    Passing these values to to lu_solve with a b vector will solve the equation"""
    dim = A.shape[0]
    pv = np.arange(0, dim, 1)
    err = 0
    # Eliminate everything but the last row (dim-1)
    for i in range(dim - 1):
        # Store the current pivot from the pivot list
        pvt = pv[i]
        # Store the value of the current pivot
        pvv = A[pvt, i]
        # Search the other row specified in the pivot list
        for k in pv:
            # Check if the other rows have larger pivot values
            val = A[k, i]
            # print("val ({0}) > pvv({1})".format(val, pvv))
            if val > pvv:
                # We found a larger row, store the value  and so we can check the others
                pvv = val
                pvt = k
        # Did we find a new pivot?
        if pvt > pv[i]:
            # If we did switch the indices in the pivot list
            # print("We switched row {0} with pivot {1} for row {2} with pivot {3}".format(pv[i], A[pv[i], i], pvt, A[pvt,i]))
            tmp = pv[i]
            pv[i] = pvt
            pv[pvt] = tmp
        # print(pv)
        # Check if the current pivot is close to 0
        # if it is, break and set the error flag
        if np.abs(A[pv[i], i]) < tol:
            err = -1
            break
        # Here we actually perform the actual elimination
        for j in range(i + 1, dim):
            # print("c = {0}/{1}".format(A[pv[j], i], A[pv[i], i]))
            c = A[pv[j], i] / A[pv[i], i]
            # print(A[pv[j], i:])
            # print((c * A[pv[i], i:]))
            A[pv[j], i:] -= (c * A[pv[i], i:])
            # print(A[pv[j], :])
            #print("Replacing index {0},{1} with value {2} with {3}".format(pv[j], i, A[pv[j], i], c))
            A[pv[j], i] = c
    # print(A)
    # Quit here if the system is singular
    if err == -1:
        return None
    else:
        return (A, pv)

def lu_solve(A, pv, b):
    """ Solves the system Ax=b given the output from lu_factor"""
    dim = A.shape[0]
    x = np.zeros(dim)
    for i in range(dim - 1):
        for j in range(i + 1, dim):
            #All of our c's are stored in A from the output of LU factor
            c = A[pv[j], i]
            #Calculate the b vector that would result from the typical elimination procedure
            b[pv[j]] -= (c * b[pv[i]])
    #print(d)
    x[-1] = b[pv[-1]] / A[pv[-1], -1]
    # Now we solve the remaining equations, this is the same as Gaussian back substitution
    # dim-2 starts means we begin at second row from the end and go until the 0th row
    for i in range(dim - 2, -1, -1):
        # Grab the corresponding b value
        sum = b[pv[i]]
        # Now we sum from the last column (dim -1 ) to the current column (i-1)
        for j in range(dim - 1, i - 1, -1):
            sum -= x[j] * A[pv[i], j]
        x[i] = sum / A[pv[i], i]
    return x
    
def inv(A, tol, err):
    """We always assume square matrices"""
    dim = A.shape[0]
    A1 = np.zeros(A.shape)
    A, pvt = lu_factor(A, tol, err)
    if err == -1:
        return None
    for i in range(dim):
        b = np.zeros(dim)
        b[i] = 1
        x = lu_solve(A, pvt, b)
        A1[:, i] = np.copy(x)
    return A1
    
def gauss_seidel(A, b, x, tol, maxi, lam):
    """ x should contain initial guesses (can be 0)"""
    dim = A.shape[0]
    #Divide everything by each row by its diagnol element
    for i in range(dim):
        tmp = A[i,i]
        for j in range(dim):
            A[i,j] /= tmp
        b[i] /= tmp
        # print(A)
    for i in range(dim):
        acc = b[i]
        for j in range(dim):
            if i == j:
                # print("Skipping i = {0} and j = {1}".format(i, j))
                continue
            else:
                acc -= A[i, j] * x[j]
        # print("Old x = {0}, new x = {1}".format(x[i], acc))
        x[i] = acc
    for i in range(maxi):
        flag = 1
        for k in range(dim):    
            acc = b[k]
            oldx = x[k]
            for j in range(dim):
                if k == j:
                    continue
                else:
                    # print('k = {0}, j={1}'.format(k, j))
                    acc -= (A[k,j] * x[j])
                    # print(acc)
            # print("Old x = {0}, new x = {1}".format(oldx, (lam * acc) + ((1-lam) * oldx)))
            x[k] = (lam * acc) + ((1-lam) * oldx)
            if flag ==1 and x[k] != 0:
                ea = abs((x[k] - oldx)/x[k]) * 100
                # print("Error is equal to {0}".format(ea))
                if ea > tol:
                    flag = 0
        if flag == 1:
            print('Breaking with ea = {0} and num iterations: {1}'.format(ea, i))
            break
    return x
    
e = 0
x2 = gauss(np.copy(A), np.copy(b), 0.001, e)
aa, pv = lu_factor(np.copy(A), 0.001, e)
x1 = lu_solve(aa, pv, np.copy(b))
print(np.dot(A,x2))
print(np.dot(A,x1))
x3 = gauss_seidel(np.copy(A), np.copy(b), np.zeros(A.shape[0]), 0.0001, 25, 1.03)
print(np.dot(A,x3))
