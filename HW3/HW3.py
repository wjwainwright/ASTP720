# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c
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

#Scipy odeint test
def pend(y, t, b, c):
    theta, omega = y
    dydt = [omega, -b*omega - c*np.sin(theta)]
    return dydt

y0 = [np.pi - 0.1, 0.0]
t = np.linspace(0, 10, 101)
sol = odeint(pend, y0, t, args=(0.25, 5.0))

plt.figure()
plt.title("Pendulum - Scipy ODEint")
plt.plot(t, sol[:, 0], 'b', label=r'$\theta(t)$')
plt.plot(t, sol[:, 1], 'g', label=r'$\omega(t)$')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.savefig("odeint.pdf")

#Comparison of my methods
def pendulum(t,y,yp,b=0.25,c=5):
    theta = y
    omega = yp
    return [omega,-b*omega - c*np.sin(theta)]

t1,theta1,omega1 = ODE.euler(pendulum,0,y0,t,full=True)
t2,theta2,omega2 = ODE.heun(pendulum,0,y0,t,full=True)
t3,theta3,omega3 = ODE.rk4(pendulum,0,y0,t,full=True)

#Euler
plt.figure()
plt.title("Pendulum - My Euler")
plt.plot(t1, theta1, 'b', label=r'$\theta(t)$')
plt.plot(t1, omega1, 'g', label=r'$\omega(t)$')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.savefig("euler.pdf")

#Heun
plt.figure()
plt.title("Pendulum - My Heun")
plt.plot(t2, theta2, 'b', label=r'$\theta(t)$')
plt.plot(t2, omega2, 'g', label=r'$\omega(t)$')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.savefig("heun.pdf")

#RK4
plt.figure()
plt.title("Pendulum - My RK4")
plt.plot(t3, theta3, 'b', label=r'$\theta(t)$')
plt.plot(t3, omega3, 'g', label=r'$\omega(t)$')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.savefig("rk4.pdf")

#Number 3
def stiff(t,y,lam=15):
    return [-lam*(y-np.cos(t))]

def stiffSol(t,lam=15):
    return lam/(1+lam**2)*( np.exp(-lam*t) + np.sin(t) + lam*np.cos(t) )

solve = np.arange(0,10.1,0.1)

#Euler
ts,ys,yps = ODE.euler(stiff,0,[1],solve,full=True)
plt.figure()
plt.title("Stiff Equation")
plt.plot(ts,ys,label='Euler',linewidth=7)
plt.plot(ts,[stiffSol(a) for a in ts],label='Solution')
plt.legend(loc='upper right')
plt.xlabel('t')
plt.grid()
plt.savefig("stiff1.pdf")


#Heun
ts,ys,yps = ODE.heun(stiff,0,[1],solve,full=True)
plt.figure()
plt.title("Stiff Equation")
plt.plot(ts,ys,label='Heun',linewidth=7)
plt.plot(ts,[stiffSol(a) for a in ts],label='Solution')
plt.legend(loc='upper right')
plt.xlabel('t')
plt.grid()
plt.savefig("stiff2.pdf")


#RK4
ts,ys,yps = ODE.rk4(stiff,0,[1],solve,full=True)
plt.figure()
plt.title("Stiff Equation")
plt.plot(ts,ys,label='RK4',linewidth=7)
plt.plot(ts,[stiffSol(a) for a in ts],label='Solution')
plt.legend(loc='upper right')
plt.xlabel('t')
plt.grid()
plt.savefig("stiff3.pdf")


#Number 4

def hydrostat(r,m,p):
    mu = 2
    rhoc = 10**5 *u.g*u.cm**-3
    G = c.G
    k = 1/20*(3/np.pi)**(2/3)*c.h**2/(c.m_e*c.u**(5/3))
    r = r*c.R_earth
    m = m*u.solMass
    p = p*u.Pa
    
    rho = mu*(p/k)**(3/5)
    
    if r == 0:
        r = 0.0001*c.R_earth
        rho = rhoc
    
    #Nightmare getting the units to work out. I hate cgs units
    return [float((4*np.pi*r**2*rho).to(u.solMass/u.earthRad)/u.solMass*u.earthRad),float((-1*G*m*rho*r**(-2)).to(u.Pa/u.earthRad)/u.Pa*u.earthRad)]

radii = np.arange(0,3.05,0.05)
#This is returning [0,0,0,0,...,0] for some reason
rad,mass,press = ODE.rk4(hydrostat,0,[0,0],radii,full=True)
#print(hydrostat(0,0,0))