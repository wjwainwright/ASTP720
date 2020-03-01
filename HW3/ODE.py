# -*- coding: utf-8 -*-

def euler(func,x0,y0,solvePoints,precision=3,full=False):
    """
    Euler method for solving ordinary differential equations
    
    Args:
        func: Takes a function y' = f(x,y)
        x0: Initial value for x
        y0: Initial value for the y array at the position x, i.e. y(x0)=y0
        solvePoints: List of x-values to return y(x) at
        precision: Order of precision, i.e. 10^(-1*precision)
        full: Boolean to return x,y,y' which is false by default
    Returns:
        Returns y(x) for all values in solvePoints by default, or the full
        arrays of x,y,y' if full=True
    """
    import numpy as np
    
    step = 10**(-1*precision)
    
    solvePoints = sorted(solvePoints)
    x = np.arange(x0,solvePoints[-1]+step,step)
    yprime = []
    y = [y0[0]]
    ans = []
    
    #Cut of 1.00000000000004 to 1.000
    x = [round(a,precision) for a in x]
    
    if not len(y0) == 1:
        yprime.append(y0[1])
    
    for i in range(len(x)-1):
        if len(y0) ==1:
            yprime.append(func(x[i],y[i])[0])
            y.append(y[i]+step*yprime[-1])
        else:
            yprime.append(yprime[i] + step*func(x[i],y[i],yprime[i])[1])
            y.append(y[i]+step*func(x[i],y[i],yprime[-1])[0])
    
    #Cut of 1.00000000000004 to 1.000
    solvePoints = [round(a,precision) for a in solvePoints]
    y = [round(a,precision) for a in y]
    yprime = [round(a,precision) for a in yprime]
    
    for val in solvePoints:
        ans.append(y[x.index(val)])
    
    if full:
        return x,y,yprime
    else:
        return ans

def heun(func,x0,y0,solvePoints,precision=3,full=False):
    """
    Heun's method for solving ordinary differential equations
    
    Args:
        func: Takes a function y' = f(x,y) or y' = f(x,y1,y2)
        x0: Initial value for x
        y0: Initial value for the y array at the position x, i.e. y(x0)=y0
        solvePoints: List of x-values to return y(x) at
        precision: Order of precision, i.e. 10^(-1*precision)
        full: Boolean to return x,y,y' which is false by default
    Returns:
        Returns y(x) for all values in solvePoints by default, or the full
        arrays of x,y,y' if full=True
    """
    
    import numpy as np
    
    step = 10**(-1*precision)
    
    solvePoints = sorted(solvePoints)
    x = np.arange(x0,solvePoints[-1]+step,step)
    yprime = []
    y = [y0[0]]
    ans = []
    
    #Cut of 1.00000000000004 to 1.000
    x = [round(a,precision) for a in x]
    
    if not len(y0) == 1:
        yprime.append(y0[1])
    
    for i in range(len(x)-1):
        
        if len(y0) == 1:
            #p=predictor
            #c=corrector
            p = y[i] + step*func(x[i],y[i])[0]
            c = y[i] + 0.5*step*(func(x[i],y[i])[0] + func(x[i+1],p)[0])
            y.append(c)
            yprime.append(func(x[i],y[i])[0])
        else:
            p = y[i] + step*func(x[i],y[i],yprime[i])[0]
            p2 = yprime[i] + step*func(x[i],y[i],yprime[i])[1]
            
            c = y[i] + 0.5*step*(func(x[i],y[i],yprime[i])[0] + func(x[i+1],p,p2)[0])
            c2 = yprime[i] + 0.5*step*(func(x[i],y[i],yprime[i])[1] + func(x[i+1],p,p2)[1])
            
            y.append(c)
            yprime.append(c2)
          
    #Cut of 1.00000000000004 to 1.000
    solvePoints = [round(a,precision) for a in solvePoints]
    y = [round(a,precision) for a in y]
    yprime = [round(a,precision) for a in yprime]
    
    for val in solvePoints:
        ans.append(y[x.index(val)])
    
    if full:
        return x,y,yprime
    else:
        return ans


def rk4(func,x0,y0,solvePoints,precision=3,full=False):
    """
    RK4 method for solving ordinary differential equations
    
    Args:
        func: Takes a function y' = f(x,y) or y' = f(x,y1,y2)
        x0: Initial value for x
        y0: Initial value for the y array at the position x, i.e. y(x0)=y0
        solvePoints: List of x-values to return y(x) at
        precision: Order of precision, i.e. 10^(-1*precision)
        full: Boolean to return x,y,y' which is false by default
    Returns:
        Returns y(x) for all values in solvePoints by default, or the full
        arrays of x,y,y' if full=True
    """
    
    import numpy as np
    
    #Step size
    step = 10**(-1*precision)
    
    solvePoints = sorted(solvePoints)
    x = np.arange(x0,solvePoints[-1]+step,step)
    yprime = []
    y = [y0[0]]
    ans = []
    
    #Cut of 1.00000000000004 to 1.000
    x = [round(a,precision) for a in x]
    
    if not len(y0) == 1:
        yprime.append(y0[1])
    
    for i in range(len(x)-1):
        
        if len(y0) == 1:
            
            k1 = func(x[i],y[i])[0]
            k2 = func(x[i]+0.5*step,y[i]+0.5*step*k1)[0]
            k3 = func(x[i]+0.5*step,y[i]+0.5*step*k2)[0]
            k4 = func(x[i+1],y[i]+step*k3)[0]
            
            y.append(y[i]+step/6*(k1+2*k2+2*k3+k4))
            yprime.append(k1)
            
        else:
            
            k1 = func(x[i],y[i],yprime[i])[0]
            k1p = func(x[i],y[i],yprime[i])[1]
            
            k2 = func(x[i]+0.5*step,y[i]+0.5*step*k1,yprime[i]+0.5*step*k1p)[0]
            k2p = func(x[i]+0.5*step,y[i]+0.5*step*k1,yprime[i]+0.5*step*k1p)[1]
            
            k3 = func(x[i]+0.5*step,y[i]+0.5*step*k2,yprime[i]+0.5*step*k2p)[0]
            k3p = func(x[i]+0.5*step,y[i]+0.5*step*k2,yprime[i]+0.5*step*k2p)[1]
            
            k4 = func(x[i+1],y[i]+step*k3,yprime[i]+step*k3p)[0]
            k4p = func(x[i+1],y[i]+step*k3,yprime[i]+step*k3p)[1]
            
            y.append(y[i]+step/6*(k1+2*k2+2*k3+k4))
            yprime.append(yprime[i]+step/6*(k1p+2*k2p+2*k3p+k4p))
    
    #Cut of 1.00000000000004 to 1.000
    solvePoints = [round(a,precision) for a in solvePoints]
    y = [round(a,precision) for a in y]
    yprime = [round(a,precision) for a in yprime]
    
    for val in solvePoints:
        ans.append(y[x.index(val)])
    
    if full:
        return x,y,yprime
    else:
        return ans
