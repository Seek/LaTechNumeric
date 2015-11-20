import numpy as np
import scipy.linalg as spla
import timeit
def gauss(a, b, allow_overwrite=False):
    """
        Solves  ax=b using gaussian elimation and returns x.
        
        Parameters:
        a: (n,n) array
        b: (n)  array
        allow_overwrite: Allow the code to overwrite the source matrices
    """
    n = a.shape[0]
    x = np.empty(n)
    pvtv = np.arange(0, n, 1)
    if allow_overwrite:
        a1 = a
        b1 = b
    else:
        a1 = np.copy(a)
        b1 = np.copy(b)
    s = np.amax(abs(a1), axis=1)
    #print('Maximum of each row is {0}'.format(s))
    #For each column
    for i in range(n-1):
        #print('Before elimination our matrix looks like')
        #print(a1)   
        gpivot(a1,b1, pvtv, i, s)
        #Check for division by 0
        if np.isclose(a[pvtv[i]/s[pvtv[i]], i], 0):
            #print('Coefficient a[{0}, {1}] with value {2} is too close to 0'.format(pvtv[i], i, a[pvtv[i], i]))
            raise ZeroDivisionError
        
        geliminate(a1, b1, pvtv, i)
        #print('After elimination our matrix looks like')
        #print(a1)
    #print('Our matrix now looks like')
    #print(a1)
    #print('Starting the substitution process')
    gsubstitute(a1, b1, pvtv, x, n-1)
    return x
    

def gpivot(a, b, pvtv, ind, s):
    #pvtv[ind] is the original row
    p = pvtv[ind]
    v = abs(a[p, ind])/s[pvtv[ind]]
    #print('The current pivot element is {0} on row {1}'.format(v, p))
    #print('Searching for a new pivot in rows {0} while excluding {1}'.format(pvtv[ind+1:], pvtv[:ind]))
    #i is the actual index in the pivot vector
    for i in range(ind+1, len(pvtv)):
        vv = abs(a[pvtv[i], ind])/s[pvtv[i]]
        #print('Comparing {0} and {1} for a new pivot'.format(vv, v))
        if vv > v:
            #print('Found a new potential pivot at row {0}'.format(pvtv[i]))
            v = vv
            p = i
    if p != pvtv[ind]:
        #print('The new pivot is row {0} and the old was {1}'.format(pvtv[p], pvtv[ind]))
        tmp = pvtv[ind]
        pvtv[ind] = pvtv[p]
        pvtv[p] = tmp
    #print('The pivot vector looks like: {0}'.format(pvtv))

def geliminate(a, b, pvtv, ind):
    for i in pvtv[ind+1:]:
        #print('Eliminating row {0} with row {1}'.format(i, pvtv[ind]))
        c = a[i, ind]/a[pvtv[ind], ind]
        a[i, ind:] -= c*a[pvtv[ind], ind:]
        b[i] -= c*b[pvtv[ind]]

def gsubstitute(a, b, pvtv, x, n):
    x[n] = b[pvtv[n]]/a[pvtv[n], n]
    #print('The first x value is {0}'.format(x[n]))
    for i in range(n-1, -1, -1):
        #print("We're on column {0} and we're starting on row {1}".format(i, pvtv[i]))
        #print('Starting with b = {0}'.format(b[i]))
        acc = b[pvtv[i]]
        for j in range(n , i, -1):
            acc -= a[pvtv[i], j] * x[j]
        x[i] = acc/a[pvtv[i], i]

if __name__ == "__main__":
    m = (-200) * np.random.rand(25,25) + 100
    mm = (-200) * np.random.rand(25) + 100
    #print("The answers are {0}".format(mm))
    x1 = spla.solve(m,mm)
    x2 = gauss(m,mm)
    print("x vector according to scipy is {0}".format(x1))
    print("x vector according to me is {0}".format(x2))