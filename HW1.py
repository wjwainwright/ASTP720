# -*- coding: utf-8 -*-

import rootFind as rf

def sq(num):
    return lambda x: x**2 - num
def sqPrime(num):
    return lambda x: 2*x

print(rf.bisect(sq(45),0,100))
print(rf.newton(sq(45),sqPrime(45),10))
print(rf.secant(sq(45),0,100))