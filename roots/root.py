import sys
import math
EPS = sys.float_info.epsilon

def f(x):
	return x**2  - 2
def df(x):
	return 2*x
	
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
	
def newton(f, df, x0, eps, maxn):
	xr = x0
	err = 100
	for i in range(maxn):
		xold = xr
		xr = xold - f(xold)/df(xold)
		if not abs(xr - 0) <  EPS:
			err = abs((xr-xold)/xr) * 100
		if err < eps:
			break
	return xr
	
newton(f, df, 1, 0.001, 100)

def secant(f, x0, x1, eps, maxn):
	xl = x0
	xr = x1
	err = 100
	for i in range(maxn):
		xro = xr
		xr = xro - (f(xro)*(xl-xro))/(f(xl)-f(xro))
		if not abs(xr - 0) <  EPS:
			err = abs((xr-xro)/xr) * 100
		if err < eps:
			break
		xl = xro
	return xr
secant(f, 0, 2, 0.001, 100)