# -*- coding: utf-8 -*-

import rootFind as rf
import interpolate as ip
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

def sq(num):
    return lambda x: x**2 - num
def sqPrime(num):
    return lambda x: 2*x


def gaussianLens(xprime,lam,re,N0,D,a):
    return lambda x: x*u.AU*( 1+ 1/(np.pi*a**2)*lam*re**2*N0*D*np.exp(-1*(x*u.AU/a)**2) ) - xprime

def pseudoIsotherm(x):
    
    return

def pseudoIsothermFWHM(x):
    return (1+x**2)**-0.5 - 0.5

def PrimeFWHM(x):
    return -1*x*(1+x**2)**-1.5



#Number 1 test
print(rf.bisect(sq(45),0,100))
print(rf.newton(sq(45),sqPrime(45),10))
print(rf.secant(sq(45),0,100))

#Number 2 test
xarr=[0,1,2,3,4,5]
yarr=[0,1,4,9,16,25]

func = ip.linearInterpolate(xarr,yarr)
print(func(3.5))

#Number 3
rt,count=rf.bisect(pseudoIsothermFWHM,0,100)
print(rt)

thresholds = 10**np.linspace(-10,-1,10)
counts_bisect = []
counts_newton = []
counts_secant = []

for thresh in thresholds:
    r,c = rf.bisect(pseudoIsothermFWHM,0,100,threshold=thresh)
    counts_bisect.append(c)
    
    r,c = rf.newton(pseudoIsothermFWHM,PrimeFWHM,1,threshold=thresh)
    counts_newton.append(c)
    
    r,c = rf.secant(pseudoIsothermFWHM,-1,10,threshold=thresh)
    counts_secant.append(c)

plt.figure("Counts")
plt.xlabel("Log10(Threshold)")
plt.ylabel("Iteration Count")
plt.plot(np.log10(thresholds),counts_bisect,label='Bisect')
plt.plot(np.log10(thresholds),counts_newton,label='Newton')
plt.plot(np.log10(thresholds),counts_secant,label='Secant')
plt.legend()
plt.savefig("Threshold_Count_Plot.png")
plt.savefig("Threshold_Count_Plot.pdf")

#Number 4

xprimelist = np.arange(0,2,0.025)*u.AU
xlist = []

lam = 21*u.cm
re = 2.817e-17*u.cm
N0 = 0.01*u.pc*(u.cm)**-3
D = 1*u.kpc
a = 1*u.AU

for xprime in xprimelist:
    func = gaussianLens(xprime,lam,re,N0,D,a)
    x,count = rf.bisect(func,0,1000,threshold=0.0001*u.AU)
    xlist.append(x)

yprimelist = np.zeros(len(xlist))
ylist = np.full(len(xprimelist),10)

plt.figure("Gaussian Lens")
plt.title("Gaussian Lens")
plt.xlabel("x'")
plt.yticks([])
for i in range(0,len(xlist)):
    xl = [xprimelist[i]/u.AU,xlist[i]]
    yl = [yprimelist[i],ylist[i]-np.exp(-1*((xlist[i]-1)/0.25)**2)] #The gaussian here is just for visual
    plt.plot(xl,yl,color='blue')