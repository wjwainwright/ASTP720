# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt


class node:
    
    def __init__(self,parent,height,xmin,xmax,ymin,ymax,zmin,zmax,letter):
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
        self.letter = letter
        
    def addChild(self,ch):
        self.child.append(ch)
        self.children += 1
        
    def addGal(self,gal):
        self.galaxies.append(gal)
        self.gals += 1
        return
    
    def galCount(self):
        
        
        count = 0
        
        #Check if node has children
        if self.children > 0:
            
            #Loop through children and recursively return their counts
            for ch in self.child:
                count += ch.galCount()
        
        #Recursive end-point which returns count in a terminal node
        else:
            count = self.gals
        return count
    
    
    def calc_COM(self):
        
        if self.children == 0:
            if self.gals > 0:
                g = self.galaxies[0]
                
                self.xcom = g.x
                self.ycom = g.y
                self.zcom = g.z
                self.mtot = g.mass
                
                return g.x,g.y,g.z,g.mass
        else:
            xcom = 0
            ycom = 0
            zcom = 0
            mtot = 0
            
            for ch in self.child:
                
                sub = ch.calc_COM()
                
                if sub == None:
                    return
                
                xcom += sub[0]*sub[3]
                ycom += sub[1]*sub[3]
                zcom += sub[2]*sub[3]
                mtot += sub[3]
                
            self.xcom = xcom / mtot
            self.ycom = ycom / mtot
            self.zcom = zcom / mtot
            self.mtot = mtot
        
            return self.xcom,self.ycom,self.zcom,self.mtot
    
    
    
class tree:
    
    def __init__(self,xmin,xmax,ymin,ymax,zmin,zmax):
        self.gals = []
        self.nodes = 0
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.maxDepth = 0
        self.root = node('Me',0,xmin,xmax,ymin,ymax,zmin,zmax,'Root')
        self.subDiv(self.root)
    
    def addGal(self,gal,pnode):
        
        if pnode == self.root:
            try:
                assert gal.x <= pnode.xmax and gal.x >= pnode.xmin and gal.y <= pnode.ymax and gal.y >= pnode.ymin and gal.z <= pnode.zmax and gal.z >= pnode.zmin
            except:
                print(f"{gal.x}   {gal.y}   {gal.z}   |   {pnode.xmin} - {pnode.xmax}   {pnode.ymin} - {pnode.ymax}   {pnode.zmin} - {pnode.zmax}")
        
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
                        for g in nd.galaxies:
                            self.addGal(g,nd)
                        nd.gals = 0
                        nd.galaxies = []
                        
                #If the current box has subdivisions already, call the method again with current child node as the new parent
                #This is what makes the loop/function recursive other than the subdivision and re-settling above
                else:
                    self.addGal(gal,nd)
    
    
    def toList(self,nd="root"):
        treeList = []
        
        if nd == "root":
            nd = self.root

        if nd.children == 0:
            return [nd]
        else:
            for ch in nd.child:
                treeList.extend(self.toList(ch))
            return treeList
    
    
    
    def plot(self):
        from matplotlib.patches import Rectangle
        
        nodeList = self.toList()
        
        plt.figure()
        plt.title("NOTE: z-axis is ignored")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.axis('square')
        
        plt.xlim(-1,11)
        plt.ylim(-1,11)
        
        for nd in nodeList:
            try:
                assert nd.xmax-nd.xmin == nd.ymax-nd.ymin
            except:
                print(nd.xmax,nd.xmin,nd.ymax,nd.ymin,nd.letter,"   |   ",nd.parent.xmin,nd.parent.xmax,nd.parent.ymin,nd.parent.ymax)
            plt.gca().add_patch(Rectangle((nd.xmin,nd.ymin),nd.xmax-nd.xmin,nd.ymax-nd.ymin,edgecolor='black',fill=None))
            
            if nd.gals > 0:
                for g in nd.galaxies:
                    plt.scatter(g.x,g.y,color='red')
            
        plt.show()
    
    
    
    
    
    def subDiv(self,parent):
        xmid = (parent.xmax+parent.xmin)/2
        ymid = (parent.ymax+parent.ymin)/2
        zmid = (parent.zmax+parent.zmin)/2
        
        x0 = parent.xmin
        x1 = parent.xmax
        y0 = parent.ymin
        y1 = parent.ymax
        z0 = parent.zmin
        z1 = parent.zmax
        
        if parent.height+1 > self.maxDepth:
            self.maxDepth = parent.height+1
        
        parent.addChild(node(parent,parent.height+1,x0,xmid,y0,ymid,z0,zmid,'a'))
        parent.addChild(node(parent,parent.height+1,xmid,x1,y0,ymid,z0,zmid,'b'))
        parent.addChild(node(parent,parent.height+1,x0,xmid,ymid,y1,z0,zmid,'c'))
        parent.addChild(node(parent,parent.height+1,xmid,x1,ymid,y1,z0,zmid,'d'))
        parent.addChild(node(parent,parent.height+1,x0,xmid,y0,ymid,zmid,z1,'e'))
        parent.addChild(node(parent,parent.height+1,xmid,x1,y0,ymid,zmid,z1,'f'))
        parent.addChild(node(parent,parent.height+1,x0,xmid,ymid,y1,zmid,z1,'g'))
        parent.addChild(node(parent,parent.height+1,xmid,x1,ymid,y1,zmid,z1,'h'))
    
    
        
        



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
cl1 = np.load('galaxies0.npy')
cl1 = np.load('galaxies1.npy')

#Box size
xmin=np.floor(sorted(cl1[:,0])[0])
xmax=np.ceil(sorted(cl1[:,0])[-1])
ymin=np.floor(sorted(cl1[:,1])[0])
ymax=np.ceil(sorted(cl1[:,1])[-1])
zmin=np.floor(sorted(cl1[:,2])[0])
zmax=np.ceil(sorted(cl1[:,2])[-1])

mass = 10**12

#Timestep
dt = 1000
time = np.arange(0,1e6,dt)

ref = tree(xmin,xmax,ymin,ymax,zmin,zmax)
series = [tree(xmin,xmax,ymin,ymax,zmin,zmax) for a in time]

for g in cl1:
    series[0].addGal(galaxy(g[0],g[1],g[2],mass),series[0].root)
print(series[0].maxDepth)
print(series[0].root.galCount())
cList = series[0].toList()
print(len(cList))
print(series[0].root.calc_COM())

series[0].plot()



"""
for t in time:
    series.append(evolve(series[-2],series[-1],dt))

"""