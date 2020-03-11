# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt


class node:
    
    def __init__(self,parent,height,xmin,xmax,ymin,ymax,zmin,zmax):
        self.parent = parent
        self.height = height
        self.child = []
        self.children = 0
        self.gals = 0
        self.galaxies = []
        
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        
    def addChild(self,child):
        self.child.append(child)
        self.children += 1
        
    def addGal(self,gal):
        self.galaxies.append(gal)
        self.gals += 1
        return
    
    def galCount(self):
        count = 0
        
        #Loop through children 
        if self.children > 0:
            for child in self.child:
                count += child.galCount()
        else:
            count = self.gals
        return count
    
class tree:
    
    def __init(self,xmin,xmax,ymin,ymax,zmin,zmax):
        self.gals = []
        self.nodes = 0
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.parent = node('Me',0,xmin,xmax,ymin,ymax,zmin,zmax)
        self.subDiv(self.parent)
    
    def addGal(self,gal,pnode):
        #Recursive loop that checks immediate children of current parent node
        for nd in pnode.child:
            
            #Check if the galaxy belongs in this box i.e this node
            if gal.x <= nd.xmax and gal.x >= nd.xmin and gal.y <= nd.ymax and gal.y >= nd.ymin and gal.z <= nd.zmax and gal.z >= nd.zmin:
                
                #Check if this box is undivided
                if nd.children == 0:
                    nd.addGal(gal)
                    
                    #If there is a second galaxy in this box, subdivide, clear the box, and re-allocate these galaxies
                    #This should recursively repeat until all galaxies are settled in their own box
                    if nd.gals > 1:
                        self.subDiv(nd)
                        for g in nd.gals:
                            self.addGal(g)
                        nd.gals = 0
                        nd.galaxies = []
                        
                #If the current box has subdivisions already, call the method again with current child node as the new parent
                #This is what makes the loop/function recursive other than the subdivision and re-settling above
                else:
                    self.addGal(gal,nd)
    
    def subDiv(self,parent):
        xmid = (parent.xmax-parent.xmin)/2
        ymid = (parent.ymax-parent.ymin)/2
        zmid = (parent.zmax-parent.zmin)/2
        
        x0 = parent.xmin
        x1 = parent.xmax
        y0 = parent.ymin
        y1 = parent.ymax
        z0 = parent.zmin
        z1 = parent.zmax
        
        
        parent.addChild(node(parent,parent.height+1,x0,xmid,y0,ymid,z0,zmid))
        parent.addChild(node(parent,parent.height+1,xmid,x1,y0,ymid,z0,zmid))
        parent.addChild(node(parent,parent.height+1,x0,xmid,ymid,y1,z0,zmid))
        parent.addChild(node(parent,parent.height+1,xmid,x1,ymid,y1,z0,zmid))
        parent.addChild(node(parent,parent.height+1,x0,xmid,y0,ymid,zmid,z1))
        parent.addChild(node(parent,parent.height+1,xmid,x1,y0,ymid,zmid,z1))
        parent.addChild(node(parent,parent.height+1,x0,xmid,ymid,y1,zmid,z1))
        parent.addChild(node(parent,parent.height+1,xmid,x1,ymid,y1,zmid,z1))
        
        



class galaxy:
    
    def __init__(self,x,y,z,mass):
        self.x = x
        self.y = y
        self.z = z
        self.mass = mass
        
"""
class cluster:
    
    def __init__(self,inp,mass,mode='data'):
        if mode == 'file':
            galaxies = []
            
            fn = np.load(inp)
            
            for g in fn:
                gal = galaxy(g[0],g[1],g[2],mass)
                galaxies.append(gal)
                
            
            
                
            self.N = len(fn)
            self.mass = mass*self.N
        else:
            self.galaxies,self.N,Mass = inp
            self.mass = mass*self.N
"""

def evolve(g0,g1,dt):
    
    return





#Create clusters
series = []

series.append(tree())



#Timestep
dt = 1000
time = np.arange(0,1e6,dt)

"""
for t in time:
    series.append(evolve(series[-2],series[-1],dt))

"""