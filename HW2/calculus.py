# -*- coding: utf-8 -*-

def differentiate(xarr,yarr):
    """
    Symmetric difference numerical differentiation function
    
    Args:
        xarr: Input array of x-values used to construct the numerical derivative
        yarr: Input array of y-values used to construct the numerical derivative
    Returns:
        Returns a funcion in the form of f(x) which in turn returns the numerical
        derivative at the point x based on the input xarr and yarr
    """
    def subFunc(x):
        try:
            i = xarr.index(x)
        except(ValueError):
            print(f"{x} is not an element of the list of points")
            return None
        
        if i == 0 or i == len(xarr) - 1 :
            print(f"Cannot perform symmetric differencing on an end point")
            return None
        
        yplus = yarr[i+1]
        yminus = yarr[i-1]
        h = xarr[i+1] - xarr[i]
        return (yplus - yminus)/(2*h)
    return subFunc

def midpointInt(func):
    """
    Midpoint method numerical integration method
    
    Args:
        func: Takes a function f(x) for the purpose of finding the numerical derivative
    Returns:
        Returns a function f(a,b) which in turn returns the numerical integral of f(x)
        over the interval (a,b) using the midpoint method
    """
    def subFunc(a,b):
        import numpy as np
        
        xarr = np.linspace(a,b,1000)
        h = abs(xarr[1]-xarr[0])
        ans = []
        for i in range(len(xarr)-1):
            ans.append(func((xarr[i]+xarr[i+1])/2)*h)
            
        return sum(ans)
    return subFunc

def trapezoidInt(func):
    """
    Trapezoid method numerical integration method
    
    Args:
        func: Takes a function f(x) for the purpose of finding the numerical derivative
    Returns:
        Returns a function f(a,b) which in turn returns the numerical integral of f(x)
        over the interval (a,b) using the trapezoid method
    """
    def subFunc(a,b):
        import numpy as np
        
        xarr = np.linspace(a,b,1000)
        yarr = [2*func(x) for x in xarr]
        yarr[0] = yarr[0]/2
        yarr[-1] = yarr[-1]/2
        h = abs(xarr[1]-xarr[0])
        
        return h/2 * sum(yarr)
    return subFunc

def simpsonInt(func):
    """
    Simpson's rule numerical integration method
    
    Args:
        func: Takes a function f(x) for the purpose of finding the numerical derivative
    Returns:
        Returns a function f(a,b) which in turn returns the numerical integral of f(x)
        over the interval (a,b) using Simpson's rule
    """
    def subFunc(a,b):
        return (b-a)/6*(func(a)+4*func((a+b)/2)+func(b))
    return subFunc


