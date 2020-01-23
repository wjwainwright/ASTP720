# -*- coding: utf-8 -*-


def bisect(func,a,b,threshold=0.0001):
    """
    Bisect root finding method
    
    Args:
        func: Input function that takes a single variable i.e. f(x) whose root you want to find
        a: lower bound of the range of your initial guess where the root is an element of [a,b]
        b: upper bound of the range of your initial guess where the root is an element of [a,b]
        threshold: degree of accuracy you want in your root such that |f(root)| < threshold
    """
    
    c = (a+b)/2
    
    while(True):
        if func(a) < 0 and func(b) < 0 :
            print("Both f(a) and f(b) are negative. Try again, buddy")
            return
        elif func(a) > 0 and func(b) > 0 :
            print("Both f(a) and f(b) are positive. Try again, buddy")
            return
        else:
            if func(a)*func(c) < 0 :
                b = float(c)
            else:
                a = float(c)
            c = (a+b)/2
            if(abs(func(c)) < threshold):
                return c


def newton(func,funcPrime,pos,threshold=0.0001):
    while(abs(func(pos)) > threshold):
        pos = pos - func(pos)/funcPrime(pos)
    return pos


def secant(func,a,b,threshold=0.0001):
    while(abs(func(b)) > threshold):
        temp = float(b)
        b = b - func(b) * ( (b-a)/(func(b)-func(a)) )
        a = float(temp)
    return b


