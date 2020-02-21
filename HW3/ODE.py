# -*- coding: utf-8 -*-

def find_nearest(array, value):
    #Imports
    import numpy as np
    
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def euler(func,x0,y0,solvePoints,precision=3,full=False):
    import numpy as np
    
    step = 10**(-1*precision)
    
    solvePoints = sorted(solvePoints)
    x = np.arange(x0,solvePoints[-1]+step,step)
    yprime = []
    y = [y0[0]]
    ans = []
    
    x = [round(a,precision) for a in x]
    
    if not len(y0) == 1:
        yprime.append(y0[1])
    
    for i in range(len(x)-1):
        if len(y0) ==1:
            yprime.append(func(x[i],y[i]))
        else:
            yprime.append(yprime[i] + step*func(x[i],y[i],yprime[i]))
        y.append(y[i]+step*yprime[-1])
    
    
    y = [round(a,precision) for a in y]
    yprime = [round(a,precision) for a in yprime]
    
    for val in solvePoints:
        ans.append(y[x.index(val)])
    
    if full:
        return x,y,yprime
    else:
        return ans

def heun(func,x0,y0,solvePoints,precision=3,full=False):
    import numpy as np
    
    step = 10**(-1*precision)
    
    solvePoints = sorted(solvePoints)
    x = np.arange(x0,solvePoints[-1]+step,step)
    yprime = []
    y = [y0[0]]
    ans = []
    
    x = [round(a,precision) for a in x]
    
    if not len(y0) == 1:
        yprime.append(y[1])
    
    for i in range(len(x)-1):
        
        if len(y0) == 1:
            p = y[i] + step*func(x[i],y[i])
            c = y[i] + 0.5*step*(func(x[i],y[i]) + func(x[i+1],p))
            y.append(c)
            yprime.append(func(x[i],y[i]))
        else:
            p = y[i] + step*func(x[i],y[i],yprime[i])
            p2 = yprime[i] + step*func(x[i],y[i],yprime[i])
            
            c = y[i] + 0.5*step*(func(x[i],y[i],yprime[i]) + func(x[i+1],p,p2)
