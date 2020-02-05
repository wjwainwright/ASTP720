# -*- coding: utf-8 -*-

import numpy as np
import calculus as calc

def testFunc(x):
    return 2*x**3 - 3*x + 7

#Number 1 test
xarr = list(np.arange(0,20,0.1))
xarr = [round(x,4) for x in xarr]
yarr = [testFunc(x) for x in xarr]

#Symmetric difference derivative
diff = calc.differentiate(xarr,yarr)
print(diff(19.9))

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