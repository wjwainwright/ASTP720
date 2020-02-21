# -*- coding: utf-8 -*-

import ODE

def fPrime(x,y1,y2):
    #print(f"{x}   {y1}   {y2}")
    return 3*x**2+8*y1-3*y2

def fPrime2(x,y):
    return 4*x+5*y+1

solvePoints = [0,1,2,3]

#Euler
print(ODE.euler(fPrime,0,[-1,0],solvePoints))
print(ODE.euler(fPrime2,0,[1],solvePoints))



