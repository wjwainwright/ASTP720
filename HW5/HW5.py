# -*- coding: utf-8 -*-

#Imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.optimize as so

#Fitting with no uncertainty
def simpFit(x,y,z):
    #Populate the design matrix
    A = np.zeros((len(x),3))
    for i in range(len(x)):
        A[i,0] = 1
        A[i,1] = x[i]
        A[i,2] = z[i]
    
    #Matrix math
    ATA = np.transpose(A) @ A  
    ATAinv = np.linalg.inv(ATA)
    ATb = np.dot(np.transpose(A),y)
    params = np.dot(ATAinv,ATb)
    
    return params,ATAinv #2nd matrix is the covariance/variance matrix

#Fitting with uncertainty/noise
def lsFit(x,y,z,N):
    #Populate the noise matrix
    noise = np.zeros((len(N),len(N)))
    for i in range(len(N)):
        noise[i,i] = N[i]
    
    #Populate the design matrix
    A = np.zeros((len(x),3))
    for i in range(len(x)):
        A[i,0] = 1
        A[i,1] = x[i]
        A[i,2] = z[i]
        
    #Matrix math is dumb, don't @ me
    ATA = (np.transpose(A) @ noise) @ A  
    ATAinv = np.linalg.inv(ATA)
    ATb = np.dot((np.transpose(A) @ noise),y)
    params = np.dot(ATAinv,ATb)
    
    return params,ATAinv #2nd matrix is the covariance/variance matrix


#Unpack the data file and cast arrays to the right types
name,period,distance,vmag,jmag,hmag,kmag,emag,feh = np.genfromtxt("cepheid_data.txt",delimiter=',',skip_header=1,unpack=True,dtype='str')
period = [float(a) for a in period]
distance = [float(a)*1000 for a in distance]
vmag = [float(a) for a in vmag]
jmag = [float(a) for a in jmag]
hmag = [float(a) for a in hmag]
kmag = [float(a) for a in kmag]
emag = [float(a) for a in emag]
feh = [float(a) for a in feh]


#Function that plots in 3D
def func(x,alpha,beta,gamma):
    logP,feh = x
    return alpha + beta*logP+gamma*feh

#Function that plots in 2D
def func2(logP,alpha,beta,gamma):
    return alpha + beta*logP


#Constants
rv = 3.1

#Math to find Mv
logP = [np.log10(a) for a in period]
av = [rv*a for a in emag]
Mv = []
for i in range(len(vmag)):
    Mv.append(vmag[i]+5-5*np.log10(distance[i])-av[i])


#Simple fit with no uncertainty
myfit,myvar = simpFit(logP,Mv,feh)


#2D plot
plt.figure()
plt.scatter(logP,Mv,label='Data')
x = np.linspace(0,2,100)
plt.plot(x,[func2(a,*myfit) for a in x],color='red',label=f'My Fit: \n$\\alpha$ = {myfit[0]:.2f}±{myvar[0,0]:.2f}\n$\\beta$ = {myfit[1]:.2f}±{myvar[1,1]:.2f}\n$\\gamma$ = {myfit[2]:.2f}±{myvar[2,2]:.2f}')
plt.title('Simple fit - No $\\sigma$')
plt.xlabel('Period (Days)')
plt.ylabel(r'$M_v$')
plt.legend(loc='upper right')
plt.savefig("simplefit.pdf")



#Full fit with noise
myfit,myvar = lsFit(logP,Mv,feh,[0.1]*len(period))


#3D plot
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.set_xlabel('Period (Days)')
ax.set_ylabel(r'$\frac{Fe}{H}$')
ax.set_zlabel(r'$M_v$')

ax.scatter(logP,feh,Mv,label='Data',color='red')

x, y = np.meshgrid(range(0,3), range(-1,2))
zs = np.array(func((np.ravel(x),np.ravel(y)),*myfit))
z = zs.reshape(x.shape)
ax.plot_surface(x,y,z,alpha=0.4,label='Scipy Fit')
plt.savefig("3Dplot.pdf")


#2D plot
plt.figure()
plt.errorbar(logP,Mv,[0.1]*len(period),fmt='o',label='Data')
x = np.linspace(0,2,100)
plt.plot(x,[func2(a,*myfit) for a in x],color='red',label=f'My Fit: \n$\\alpha$ = {myfit[0]:.2f}±{myvar[0,0]:.2f}\n$\\beta$ = {myfit[1]:.2f}±{myvar[1,1]:.2f}\n$\\gamma$ = {myfit[2]:.2f}±{myvar[2,2]:.2f}')
plt.title('Full Fit - $\\sigma$=0.1')
plt.xlabel('Period (Days)')
plt.ylabel(r'$M_v$')
plt.legend(loc='upper right')
plt.savefig("fullfit.pdf")

#Scipy optimize fit
fit,var = so.curve_fit(func,(logP,feh),Mv)


#Scipy comparison
plt.figure()
plt.errorbar(logP,Mv,[0.1]*len(period),fmt='o',label='Data')
x = np.linspace(0,2,100)
plt.plot(x,[func2(a,*myfit) for a in x],color='red',label=f'My Fit: \n$\\alpha$ = {myfit[0]:.2f}±{myvar[0,0]:.2f}\n$\\beta$ = {myfit[1]:.2f}±{myvar[1,1]:.2f}\n$\\gamma$ = {myfit[2]:.2f}±{myvar[2,2]:.2f}')
plt.plot(x,[func2(a,*fit) for a in x],'--g',label=f'Scipy Fit: \n$\\alpha$ = {fit[0]:.2f}±{var[0,0]:.2f}\n$\\beta$ = {fit[1]:.2f}±{var[1,1]:.2f}\n$\\gamma$ = {fit[2]:.2f}±{var[2,2]:.2f}')
plt.title('Scipy Comparison - $\\sigma$=0.1')
plt.xlabel('Period (Days)')
plt.ylabel(r'$M_v$')
plt.legend(loc='upper right')
plt.savefig("comparison.pdf")


