# -*- coding: utf-8 -*-

def linearInterpolate(xarr,yarr):
    """
    Linear piecewise interpolation function
    
    Args:
        xarr: Array of x values that is in the same order as the y value array
        yarr: Array of y values that is in the same order as the x value array
        
    Returns:
        Returns a function in the form f(x) which, when given an x value inside
        the initial bounds of xarr and yarr, will return a y value according to
        a linear piecewise interpolation. If x is out of bounds it returns none
    """
    
    def subFunc(x):
        import numpy as np
        
        if x < xarr[0] or x > xarr[-1]:
            print("X position out of bounds")
            return
        
        arr = np.asarray(xarr) 
        a = arr[(np.abs(arr - x)).argmin()]
        
        if a < x:
            x0 = a
            x1 = xarr[xarr.index(a)+1]
        else:
            x1 = a
            x0 = xarr[xarr.index(a)-1]
        y0 = yarr[xarr.index(x0)]
        y1 = yarr[xarr.index(x1)]
        xd = (x-x0)/(x1-x0)
        
        return y0*(1-xd)+y1*xd
        
    return subFunc