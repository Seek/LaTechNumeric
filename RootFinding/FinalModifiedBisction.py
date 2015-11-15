import sys
EPS = sys.float_info.epsilon


#Define the function
def f(x):
	return (x+1)**2 - 1
	
def bisect(f, x1, x2, eps, maxn):
	assert f(x1)*f(x2) < 0, \
	"We cannot find a root if the function does not change signs"
	xl = x1
	xu = x2
	xr = 0
	fl = f(xl)
	err = 1000
	for i in range(maxn):
		r = (xl + xu)/2
		print(r)
		fr = f(r)
		if not abs(r - 0) <  EPS:
			err = abs((r-xr)/r) * 100
		if err < eps:
			print("Error =" + str(err))
			break
		v = fl * fr
		if v < 0:
			xu = r
		if v > 0:
			xl = r
			fl = fr
		else:
			err = 0
		xr = r
	return r
	
print("Computing the roots of x**2 - 2")
r = bisect(f, -1.5, 10, 0.00001, 100)
print("Root = " + str(r))

