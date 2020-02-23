# -*- coding: utf-8 -*-

import ODE

def f2Prime(x,y1,y2):
    #print(f"{x}   {y1}   {y2}")
    return [3*x**2+8*y1-3*y2]*2

def fPrime(x,y):
    return [5*y+1]

solvePoints = [0,1,2,3]

#Euler
print(ODE.euler(f2Prime,0,[-1,0],solvePoints))
print(ODE.euler(fPrime,0,[-1],solvePoints))

#Heun
print(ODE.heun(f2Prime,0,[-1,0],solvePoints))
print(ODE.heun(fPrime,0,[-1],solvePoints))

#RK4
print(ODE.rk4(f2Prime,0,[-1,0],solvePoints))
print(ODE.rk4(fPrime,0,[-1],solvePoints))
