# -*- coding: utf-8 -*-

import numpy as np

class matrix:
    
    def __init__(self,m,n,source='none',mode='none'):
        self.m = int(m)
        self.n = int(n)
        self.values = np.zeros((self.m,self.n))
        
        if not mode == 'none':
             
            if mode == 'numpy':
                #Importing a numpy matrix to my matrix format (which is still a numpy matrix under the hood)
                self.values = source
            else:
                file = open(source).read()
                lines = file.split('\n')
                lines.pop(0)
                
                for line in lines:
                    line = line.strip()
                    val = line.split(',')
                    
                    r = int(val[0]) - 1
                    c = int(val[1]) - 1
                    self.values[r,c] = float(val[2])
                
    def __add__(self,other):
        
        try:
            assert self.m == other.m
            assert self.n == other.m
        except:
            print("Matrices must be the same shape to add")
            raise
        
        
        ans = matrix(self.m,self.n)
        
        for r,row in enumerate(self.values):
            for c,val in enumerate(row):
                ans.values[r,c] = val + other.values[r,c]
        return ans
    
    def __sub__(self,other):
        
        try:
            assert self.m == other.m
            assert self.n == other.m
        except:
            print("Matrices must be the same shape to subtract")
            raise
        
        
        ans = matrix(self.m,self.n)
        for r,row in enumerate(self.values):
            for c,val in enumerate(row):
                ans.values[r,c] = val - other.values[r,c]
        return ans
    
    def __mul__(self,other):
        try:
            assert self.n == other.m
        except:
            print("Matrices must have compatible dimensions for multiplication. M1_cols == M2_rows")
            raise
        
        ans = matrix(self.m,other.n)
        
        for r,row in enumerate(self.values):
            for c,col in enumerate(other.transpose().values):
                ans.values[r][c] = sum([a*b for a,b in zip(row,col)])
                
        return ans

    def transpose(self):
        ans = matrix(self.n,self.m)
        
        for r,row in enumerate(self.values):
            for c,val in enumerate(row):
                ans.values[c,r] = val
        return ans
    
    def trace(self):
        try:
            assert self.m == self.n
        except:
            print("Matrix must be square")
            raise
        
        return sum([self.values[i,i] for i in range(self.m)])
    
    def determinant(self):
        try:
            assert self.m == self.n
        except:
            print("Matrix must be square")
            raise
        
        if self.m == 2:
            return self.values[0,0]*self.values[1,1] - self.values[0,1]*self.values[1,0]
        
        det = []
        for c,val in enumerate(self.values[0]):
            det.append( (-1)**c * val * self.residual(0,c).determinant() )
        return sum(det)

    def residual(self,r,c):
        ans = np.delete(self.values,r,0)
        ans = np.delete(ans,c,1)
        return matrix(self.m-1,self.n-1,ans,mode='numpy')
    
    def invert(self):
        ans = matrix(self.m,self.n)
        det = self.determinant()
        
        for r in range(self.m):
            for c in range(self.n):
                ans.values[r,c] = (-1)**(r+c) * self.residual(r,c).determinant() / det
        
        return ans
    
    def luDecomp(self):
        l = matrix(self.m,self.n)
        u = matrix(self.m,self.n)
        
        for j in range(self.n):
            
            l.values[j,j] = 1
            
            for i in range(j+1):
                u.values[i,j] = 1/l.values[i,i] * ( self.values[i,j] - sum([l.values[i,k]*u.values[k,j] for k in range(i)]) )
            for i in range(j,self.n):
                l.values[i,j] = 1/u.values[j,j] * ( self.values[i,j] - sum([l.values[i,k]*u.values[k,j] for k in range(j)]) )
            
        return l,u

    def luDeterminant(self):
        l,u = self.luDecomp()

        return np.prod([u.values[i,i] for i in range(u.n)])
    
        
    