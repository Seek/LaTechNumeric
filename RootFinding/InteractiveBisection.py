import sys
from sympy import sympify, lambdify
from sympy.abc import x
EPS = sys.float_info.epsilon
	
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

eq = input("Please enter a function of x: ")
eq = sympify(eq)
f = lambdify(x, eq, 'math')
ll = float(input("Please enter the lower limit: "))
ul = float(input("Please enter the upper limit: "))
print("Calculating...")
root = bisect(f, ll, ul, 0.001, 100)
print("The root in the given range seems to be:" + str(root))