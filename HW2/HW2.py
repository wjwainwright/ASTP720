# -*- coding: utf-8 -*-

import numpy as np
import calculus as calc
from matrix import matrix

def testFunc(x):
    return 2*x**3 - 3*x + 7

#Number 1 test
print("Number 1 Test:")
xarr = list(np.arange(0,20,0.1))
xarr = [round(x,4) for x in xarr]
yarr = [testFunc(x) for x in xarr]

#Symmetric difference derivative
diff = calc.differentiate(xarr,yarr)
print(diff(19.5))

a=-4
b=-20

#Midpoint Integral
mp = calc.midpointInt(testFunc)
print(mp(a,b))

#Trapezoid Integral
trap = calc.trapezoidInt(testFunc)
print(trap(a,b))

#Simpson Integral
sim = calc.simpsonInt(testFunc)
print(sim(a,b))





#Number 3 test
print("Number 3 Test:")
#m = matrix(9,9,source='A_coefficients.dat',mode='file')
m = np.empty((0,4))
m = np.r_[m,[[2,10,6,9]]]
m = np.r_[m,[[3,-1,4,20]]]
m = np.r_[m,[[5,0,3,-11]]]
m = np.r_[m,[[-7,8,0,2]]]

m = matrix(4,4,m,mode='numpy')

print(m.values)
summed = m+m
print(summed.values)
transpose = m.transpose()
print(transpose.values)
mult = m*m
print(mult.values)
print(m.trace())

print(m.determinant())