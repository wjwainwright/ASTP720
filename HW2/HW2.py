# -*- coding: utf-8 -*-

import numpy as np
import calculus as calc
from matrix import matrix
import astropy.units as units
import astropy.constants as const
import matplotlib.pyplot as plt

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

#2

def vc(c,v200,r200):
    #Vc(r) method as defined in eq. 2
    def subFunc(r):
        x = r/r200
        
        num = np.log(1+c*x) - (c*x)/(1+c*x)
        div = np.log(1+c) - c/(1+c)
        
        return v200*np.sqrt(num/div)*x**(-0.5)
        
    return subFunc

def menc(G,func):
    def subFunc(r):
        return r*func(r)**2/G
    return subFunc

def rho(rs):
    def subFunc(r):
        return 1/( (r/rs)*(1+r/rs)**2 )
    return subFunc
    

#c = 15
#v200 = 160 * units.km / units.s
#rs = r200/c

cList = [15,20]
v200List = [150,200] * units.km / units.s

r200 = 230 * units.kpc
G = 6.67e-11 * units.m**3 * units.kg**-1 * units.s**-2

for c in cList:
    for v200 in v200List:
        rs = r200/c
        
        rList = np.arange(1,350,1) * units.kpc
        func_vc = vc(c,v200,r200)
        vcList = [func_vc(r) for r in rList]
        func_rho = rho(rs)
        rhoList = [func_rho(r) for r in rList]
        func_menc = menc(G,func_vc)
        mencList = [func_menc(r) for r in rList]
        
        pointList = np.arange(10,350,10) * units.kpc
        shell = 1 * units.kpc
        MrList = [calc.simpsonInt(func_menc)(point-shell,point+shell) for point in pointList]
        
        #Vc
        plt.figure()
        plt.title(f"c={c}   v200={v200}   r200={r200}")
        plt.xlabel("r (kpc)")
        plt.ylabel("Vc (km/s)")
        plt.plot([r/units.kpc for r in rList],[vc/units.km*units.s for vc in vcList])
        plt.savefig(f"Vc_r_{c}_{int(v200/units.km*units.s)}.pdf")
        
        
        #Menc
        plt.figure()
        plt.title(f"c={c}   v200={v200}   r200={r200}")
        plt.xlabel("r (kpc)")
        plt.ylabel("Menc (Solar Masses)")
        plt.plot([r/units.kpc for r in rList],[menc/units.solMass for menc in mencList])
        plt.savefig(f"Menc_r_{c}_{int(v200/units.km*units.s)}.pdf")
        
        #M(r+delta_r)
        plt.figure()
        plt.title(f"c={c}   v200={v200}   r200={r200}")
        plt.xlabel("Shell Radius +/- 1 (kpc)")
        plt.ylabel("Menc (Solar Masses)")
        plt.plot([point/units.kpc for point in pointList],[mr/units.solMass/units.kpc for mr in MrList])
        plt.savefig(f"Mass_shell_{c}_{int(v200/units.km*units.s)}.pdf")
        
        #M_total
        
        
    

#Number 3 test
print("Number 3 Test:")

m = np.empty((0,4))
m = np.r_[m,[[2,10,6,9]]]
m = np.r_[m,[[3,-1,4,20]]]
m = np.r_[m,[[5,1,3,-11]]]
m = np.r_[m,[[-7,8,1,2]]]

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

l,u = m.luDecomp()

print(l.values)
print(u.values)
print(u.determinant())
print(m.luDeterminant())
assert m.luDeterminant() == m.determinant()
inv = m.invert()
print(inv.values)

#5
print("Number 5:")
A_ul = matrix(9,9,source='A_coefficients.dat',mode='file')
c = const.c 
h = const.h
nu = 0.5e15 * units.hertz
B_ul = matrix(9,9,source=A_ul.values * 2*h*nu**3 / (c**2),mode='numpy')
print(B_ul.values)