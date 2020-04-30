# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

plt.close('all')

hjd,flux = np.genfromtxt('lightcurve_data.txt',skip_header=18,unpack=True)

plt.figure()
plt.scatter(hjd,flux,s=1)
plt.ylim(0.992,1.001)
plt.xlabel('Heliocentric Julian Date (Days)')
plt.ylabel('Flux')
plt.title('Raw Data')
plt.savefig('raw.pdf')




thresh = 0.999

indices = []

#Identify the local minima
i = 0
while i < len(flux):
    if flux[i] < thresh:
        j = i
        cut = []
        while flux[j] < thresh:
            cut.append(flux[j])
            j += 1
            if j >= len(flux):
                break
        i = j
        minimum = sorted(cut)[0]
        indices.append(list(flux).index(minimum))
    else:
        i += 1





period = []
def nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


for i in range(len(flux)):
    
    idx = nearest(indices,i)
    
    #Phase scaled to hours
    period.append(24*(hjd[i]-hjd[idx]))


#Sort ordered pairs by x position
linked = zip(period,flux)
linked = sorted(linked,key=lambda x: x[0])
period,flux = zip(*linked)


plt.figure()
plt.scatter(period,flux,s=5)
plt.xlim(-7.5,7.5)
plt.ylim(0.992,1.001)
plt.xlabel('Phase (Hours)')
plt.ylabel('Flux')
plt.title('Folded and Phase Adjusted')
plt.savefig('folded.pdf')



mins = [hjd[a] for a in indices]
gaps = []

for i in range(1,len(mins)):
    gaps.append(mins[i]-mins[i-1])
AvgPeriod = np.mean(gaps)



def box(x,width,depth,center):
    if x < center + width and x > center - width:
        return 1 - depth
    else:
        return 1


plt.plot(period,[box(x,2,0.007,0) for x in period],color='red',label='Example Plot\nWidth=2\nDepth=0.007\nCenter=0')
plt.legend()
plt.savefig('test.pdf')

def mcmc(data,init,samples=1000,proposal_width=0.001):
    
    posterior = np.empty((0,len(init)))
    period,flux = data
    
    var_current = []
    for var in init:
        var_current.append(var)
    
    for sample in range(samples):
        
        post = []
        
        for i in range(len(init)):
            
            
            var_proposed = norm(var_current[i],proposal_width).rvs()
            
            params_prop = [a for a in var_current]
            params_prop[i] = var_proposed
            
            
            #Assuming sigma=1 for the simplest case
            box_current = np.prod([box(a,*var_current) for a in period])
            box_proposed = np.prod([box(a,*params_prop) for a in period])
            
            
            likelihood_current = norm(box_current,1).pdf(flux)
            likelihood_proposed = norm(box_proposed,1).pdf(flux)
            
            p_current = np.sum([np.log(a) for a in likelihood_current])
            p_proposed = np.sum([np.log(a) for a in likelihood_proposed])
            #print(p_proposed/p_current)
            accept = p_proposed / p_current > np.random.rand()
            
            if accept:
                #print(f"{i}  |  {var_current}  ->  {var_proposed}")
                var_current[i] = var_proposed
            post.append(var_current[i])
        posterior = np.r_[posterior,[post]]
    return posterior


    
post = mcmc([period,flux],[2,0.007,0],samples=1000)
fitWidth = np.mean(post[:,0])
fitDepth = np.mean(post[:,1])
fitCenter = np.mean(post[:,2])

plt.figure()
plt.scatter(period,flux,s=5)
plt.xlim(-7.5,7.5)
plt.ylim(0.992,1.001)
plt.xlabel('Phase (Hours)')
plt.ylabel('Flux')
plt.title('Best Fit parameters')
plt.plot(period,[box(x,fitWidth,fitDepth,fitCenter) for x in period],color='red',label=f'Fit Parameters\nWidth={fitWidth:.4f}\nDepth={fitDepth:.4f}\nCenter={fitCenter:.4f}')
plt.legend()
plt.savefig('fit.pdf')


plt.figure()
plt.title('Width Histogram')
plt.hist(post[:,0])
plt.savefig('histw.pdf')
plt.figure()
plt.title('Depth Histogram')
plt.hist(post[:,1])
plt.savefig('histd.pdf')
plt.figure()
plt.title('Center Histogram')
plt.hist(post[:,2])
plt.savefig('histc.pdf')




