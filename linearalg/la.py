import numpy as np
import scipy.linalg as spla

def gen_diag_dom(n, a, b):
    """ Generates a diagonally dominant random matrix"""
    m = (b-a) * np.random.rand(n,n) + a
    for i in range(n):
        tmp = np.sum(np.abs(m[i,:]))
        m[i,i] = (tmp-(10*tmp)) * np.random.rand() + tmp
    return m
def gen_tridiag(n, a, b):
    m = np.zeros((n, n))
    m[0, 0:2] = (b-a) * np.random.rand(2) + a
    m[n-1, n-2:n] = (b-a) * np.random.rand(2) + a
    for i in range(1, n-1):
        m[i, i-1:i+2] = (b-a)* np.random.rand(3) + a
    return m
    
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
        gpivot(a1, pvtv, i, s)
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
    

def gpivot(a, pvtv, ind, s):
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

def lueliminate(a, pvtv, ind):
    for i in pvtv[ind+1:]:
        #print('Eliminating row {0} with row {1}'.format(i, pvtv[ind]))
        c = a[i, ind]/a[pvtv[ind], ind]
        a[i, ind:] -= c*a[pvtv[ind], ind:]
        # print('Putting the value of c ({3}) at index a[{0},{1}] the current value is {2}'.format(i, ind, a[i,ind], c))
        a[i, ind] = c
        # print('Normally we would perform b[{0}] -= c*b[{1}]'.format(i, pvtv[ind]))

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
def lu_factor(a, allow_overwrite=False):
    """
        Factors a such that the returned matrix and pvtv table can be passed to lu_solve to solve for a given
        b vector.
        
        Parameters:
        a: (n,n) array
        allow_overwrite: Allow the code to overwrite the source matrices
    """
    n = a.shape[0]
    pvtv = np.arange(0, n, 1)
    if allow_overwrite:
        a1 = a
    else:
        a1 = np.copy(a)
    s = np.amax(abs(a1), axis=1)
    #print('Maximum of each row is {0}'.format(s))
    #For each column
    for i in range(n-1):
        #print('Before elimination our matrix looks like')
        #print(a1)   
        gpivot(a1, pvtv, i, s)
        #Check for division by 0
        if np.isclose(a[pvtv[i]/s[pvtv[i]], i], 0):
            #print('Coefficient a[{0}, {1}] with value {2} is too close to 0'.format(pvtv[i], i, a[pvtv[i], i]))
            raise ZeroDivisionError
        
        lueliminate(a1, pvtv, i)
        #print('After elimination our matrix looks like')
        #print(a1)
    #print('Our matrix now looks like')
    #print(a1)
    #print('Starting the substitution process')
    return (a1,pvtv)
def lu_solve(a, b, pvtv, allow_overwrite=False):
        n = a.shape[0]
        x = np.empty(n)
        if allow_overwrite:
            a1 = a
            b1 = b
        else:
            a1 = np.copy(a)
            b1 = np.copy(b)
        #Eliminate the b's as we would have before
        for i in range(n-1):
            #i is our current column
            for j in range(i+1, n):
                # print('Grabbing c value from a[{0},{1}]'.format(pvtv[j], i))
                # print('Performing b[{0}] - c*b[{1}] with c = {2}]'.format(pvtv[j], pvtv[i], a[pvtv[j], i]))
                b1[pvtv[j]] -= a1[pvtv[j], i] * b1[pvtv[i]]
        # print('The b1 vector before after forward looks like {0}'.format(b1))
        #Now back substitute
        gsubstitute(a1, b1, pvtv, x, n-1)
        return x
def inv(a):
    """ Returns the inverse of the matrix a"""
    aa, pvt = lu_factor(a)
    n = a.shape[0]
    ai = np.empty(a.shape)
    for i in range(n):
        b = np.zeros(n)
        b[i] = 1
        ai[:,i] = lu_solve(aa, b, pvt)
    return ai
    
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
def thomas_factor(a, allow_overwrite=False):
    """ Returns a LU matrix for a tridiagnol system"""
    n = a.shape[0]
    if allow_overwrite:
        a1 = a
    else:
        a1 = np.copy(a)
    for i in range(1,n):
        a1[i,i-1] /= a1[i-1, i-1]
        a1[i,i] -= a1[i,i-1] * a1[i-1, i]
    return a1
def thomas_solve(a, b, allow_overwrite=False):
    """Solves the factored matrix from thomas_factor for  given b vector"""
    n = a.shape[0]
    x = np.empty(n)
    if allow_overwrite:
        a1 = a
        b1 = b
    else:
        a1 = np.copy(a)
        b1 = np.copy(b)
    #Forward substitute
    for i in range(1,n):
        b1[i] -= a1[i, i-1] * b1[i-1]
    x[n-1] = b1[n-1]/a1[n-1, n-1]
    for i in range(n-2, -1, -1):
        x[i] = (b1[i]-a[i, i+1] * x[i+1])/ a1[i, i]
    return x
if __name__ == "__main__":
    n = 40
    # m = (-200) * np.random.rand(3,3) + 100
    # mm = (-200) * np.random.rand(3) + 100
    # #print("The answers are {0}".format(mm))
    # x1 = spla.solve(m,mm)
    # x2 = gauss(m,mm)
    # aa, pvt = lu_factor(m)
    # x3 = lu_solve(aa, mm, pvt)
    # print("x vector according to scipy is {0}".format(x1))
    # print("x vector according to my gauss is {0}".format(x2))
    # print('x vector according my lu is {0}'.format(x3))
    # print('Scipy gauss is close to my gauss? {0}'.format(np.isclose(x1,x2)))
    # print('Scipy gauss is close to my LU? {0}'.format(np.isclose(x1,x3)))
    a1 = gen_diag_dom(n,-100,100)
    b1 = (-200) * np.random.rand(n) + 100
    # ai = inv(m)
    # print('A * A^-1 = \\n{0}'.format(np.dot(m, ai)))
    x1 = spla.solve(a1,b1)
    x2 = gauss(a1,b1)
    aa, pvt = lu_factor(a1)
    x3 = lu_solve(aa, b1, pvt)
    x4 = gauss_seidel(a1,b1, np.zeros(n), 0.0001, 25, 1)
    # print("x vector according to scipy is {0}".format(x1))
    # print("x vector according to my gauss is {0}".format(x2))
    # print('x vector according my lu is {0}'.format(x3))
    # print('x vector according to my GS is {0}'.format(x4))
    a2 = gen_tridiag(n, -100, 100)
    print(a2)
    aa2 = thomas_factor(a2)
    print(aa2)
    x5 = thomas_solve(aa2, b1)
    x6 = spla.solve(a2,b1)
    print("x vector according to scipy is {0}".format(x6))
    print('The results of the thomas algorithm are: {0}'.format(x5))