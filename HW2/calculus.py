# -*- coding: utf-8 -*-

def differentiate(xarr,yarr):
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
    def subFunc(a,b):
        return
    return subFunc


