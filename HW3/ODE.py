# -*- coding: utf-8 -*-

def find_nearest(array, value):
    #Imports
    import numpy as np
    
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def euler(func,x0,y0,solvePoints,step=0.05):
    import numpy as np
    
    solvePoints = sorted(solvePoints)
    x = np.arange(solvePoints[0],solvePoints[-1]+step,step)
    
    
    
    
    
    return


