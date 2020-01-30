# -*- coding: utf-8 -*-


def bisect(func,a,b,threshold=0.0001):
    """
    Bisect root finding method
    
    Args:
        func: Input function that takes a single variable i.e. f(x) whose root you want to find
        a: lower bound of the range of your initial guess where the root is an element of [a,b]
        b: upper bound of the range of your initial guess where the root is an element of [a,b]
        threshold: degree of accuracy you want in your root such that |f(root)| < threshold
        
    Returns:
            If an even number (including zero) of roots is located between the points a and b, then the function returns nothing and prints accordingly.
            Otherwise, the function returns a float c where f(c) is within the threshold of zero.
    """
    
    c = (a+b)/2
    count = 0
    
    while(True):
        if func(a) < 0 and func(b) < 0 :
            print("Both f(a) and f(b) are negative. Try again, buddy")
            return 0,0
        elif func(a) > 0 and func(b) > 0 :
            print("Both f(a) and f(b) are positive. Try again, buddy")
            return 0,0
        else:
            if func(a)*func(c) < 0 :
                b = float(c)
            else:
                a = float(c)
            c = (a+b)/2
            count += 1
            if(abs(func(c)) < threshold):
                return c,count


def newton(func,funcPrime,pos,threshold=0.0001):
    """
    Newton root finding method
    
    Args:
        func: Input function that takes a single variable i.e. f(x) whose root you want to find
        funcPrime: Input function for the analytical derivative of the same function f(x)
        pos: initial guess for an x value that is ideally somewhat close to the root
        threshold: degree of accuracy you want in your root such that |f(root)| < threshold
        
    Returns:
            Returns a float c where f(c) is within the threshold of zero.
    """
    
    count = 0
    while(abs(func(pos)) > threshold):
        pos = pos - func(pos)/funcPrime(pos)
        count += 1
    return pos,count


def secant(func,a,b,threshold=0.0001):
    """
    Secant root finding method
    
    Args:
        func: Input function that takes a single variable i.e. f(x) whose root you want to find
        a: lower bound of the range of your initial guess where the root is an element of [a,b]
        b: upper bound of the range of your initial guess where the root is an element of [a,b]
        threshold: degree of accuracy you want in your root such that |f(root)| < threshold
        
    Returns:
            Returns a float c where f(c) is within the threshold of zero.
    """
    
    count = 0
    while(abs(func(b)) > threshold):
        temp = float(b)
        b = b - func(b) * ( (b-a)/(func(b)-func(a)) )
        a = float(temp)
        count += 1
    return b,count


