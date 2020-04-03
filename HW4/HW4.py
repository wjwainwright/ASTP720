# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt


class node:
    
    def __init__(self,parent,height,xmin,xmax,ymin,ymax,zmin,zmax,letter):
        """
        Initialization function for a node. Each node has a link to its parent and potential children.
        
        Args:
            parent: expects a pointer to a parent node object
            height: integer diagnostic variable that reports the depth of the node from the tree root
            xmin: minimum x value for the box-shaped node
            xmax: maximum x value for the box-shaped node
            ymin: minimum y value for the box-shaped node
            ymax: maximum y value for the box-shaped node
            zmin: minimum z value for the box-shaped node
            zmax: maximum z value for the box-shaped node
            letter: diagnostic variable for tracing issues that may arise during the node creation process
        """
        
        self.parent = parent
        self.height = height
        self.child = []
        self.children = 0
        self.gals = 0
        self.galaxies = []
        
        self.xcom = None
        self.ycom = None
        self.zcom = None
        
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.letter = letter
        
    def addChild(self,ch):
        """
        Method for adding a child to the node
        
        Args:
            ch: expects a pointer to a child node object that has been initialized
        """
        
        self.child.append(ch)
        self.children += 1
        
    def addGal(self,gal):
        """
        Method for adding a galaxy to the node
        
        Args:
            gal: expects a pointer to a galaxy object
        """
        
        self.galaxies.append(gal)
        self.gals += 1
        return
    
    
    def contains(self,gal):
        """
        Method for checking if a galaxy is contained within the boundaries of a node,
        but not necessarily contained in the node's list of galaxies yet
        
        Args:
            gal: expects a pointer to a galaxy object
        """
        
        if gal.x <= self.xmax and gal.x >= self.xmin and gal.y <= self.ymax and gal.y >= self.ymin and gal.z <= self.zmax and gal.z >= self.zmin:
            return True
        return False
    
    
    def galCount(self):
        """
        Method for recusrively counting the number of galaxies in a node and
        its tree of children nodes, returning the total count
        
        Returns:
            returns the total node count in all child nodes
        """
        
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
        """
        Method for calculating the center of mass for a node and its children
        
        Returns:
            returns a list of centers of masses to be summed up and assigns
            the center of mass to each node as it passes all the way up
        """
        
        #Leaf node
        if self.children == 0:
            
            #If the leaf node has a galaxy, center of mass is the galaxy
            if self.gals > 0:
                g = self.galaxies[0]
                
                self.xcom = g.x
                self.ycom = g.y
                self.zcom = g.z
                self.mtot = g.mass
                
                return [g.x],[g.y],[g.z],[g.mass]
            
            #Empty box
            else:
                return [],[],[],[]
        
        #Not a leaf node
        else:
            xcom = []
            ycom = []
            zcom = []
            mtot = []
            
            #Recursively call this method for the children and sum them up
            for ch in self.child:
                
                sub = ch.calc_COM()
                
                xcom.extend(a*b for a,b in zip(sub[0],sub[3]))
                ycom.extend(a*b for a,b in zip(sub[1],sub[3]))
                zcom.extend(a*b for a,b in zip(sub[2],sub[3]))
                mtot.extend(sub[3])
                
            self.xcom = sum(xcom) / sum(mtot)
            self.ycom = sum(ycom) / sum(mtot)
            self.zcom = sum(zcom) / sum(mtot)
            self.mtot = sum(mtot)
        
            return [self.xcom],[self.ycom],[self.zcom],[self.mtot]
    
    
    
class tree:
    
    def __init__(self,xmin,xmax,ymin,ymax,zmin,zmax):
        """
        Initialization method for a tree that is made up of node objects
        
        Args:
            xmin: minimum x value for the box-shaped root node
            xmax: maximum x value for the box-shaped root node
            ymin: minimum y value for the box-shaped root node
            ymax: maximum y value for the box-shaped root node
            zmin: minimum z value for the box-shaped root node
            zmax: maximum z value for the box-shaped root node
            
        """
        
        self.gals = []
        self.nodes = 0
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.galaxies = []
        self.maxDepth = 0
        #Root node is the start of the node tree
        self.root = node('Me',0,xmin,xmax,ymin,ymax,zmin,zmax,'Root')
        self.subDiv(self.root)
    
    def addGal(self,gal,pnode):
        """
        Method to add galaxies to a node/tree recursively
        """
        
        if pnode == self.root:
            
            self.galaxies.append(gal)
            
            #If the galaxy falls outside the bounding box of the root node, we have a problem
            try:
                assert gal.x <= pnode.xmax and gal.x >= pnode.xmin and gal.y <= pnode.ymax and gal.y >= pnode.ymin and gal.z <= pnode.zmax and gal.z >= pnode.zmin
            except:
                print(f"{gal.x}   {gal.y}   {gal.z}   |   {pnode.xmin} - {pnode.xmax}   {pnode.ymin} - {pnode.ymax}   {pnode.zmin} - {pnode.zmax}")
        
        #Recursive loop that checks immediate children of current parent node
        for nd in pnode.child:
            
            #Check if the galaxy belongs in this box i.e this node
            if nd.contains(gal):
                
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
        """
        Method for recursively creating a list from the tree structure
        
        Args:
            nd: expects a pointer to a node object that will act as the base
            for the tree, where the tree's root is the default
        Returns:
            returns a list containing every node below nd in the tree but not nd
        """
        
        treeList = []
        
        if nd == "root":
            nd = self.root

        if nd.children == 0:
            return [nd]
        else:
            for ch in nd.child:
                treeList.extend(self.toList(ch))
            return treeList
    
    
    
    def plotXY(self):
        """
        Plots the XY space for the entire tree, outlining nodes/boxes in black and
        plotting galaxies as red dots. Condenses the z-axis so red dots may share a
        box in the collapsed XY view but do not in the XYZ space as expected.
        """
        
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
            #Debug information for when I had boxes with non-square dimensions
            try:
                assert nd.xmax-nd.xmin == nd.ymax-nd.ymin
            except:
                print(nd.xmax,nd.xmin,nd.ymax,nd.ymin,nd.letter,"   |   ",nd.parent.xmin,nd.parent.xmax,nd.parent.ymin,nd.parent.ymax)
            plt.gca().add_patch(Rectangle((nd.xmin,nd.ymin),nd.xmax-nd.xmin,nd.ymax-nd.ymin,edgecolor='black',fill=None))
            
            if nd.gals > 0:
                for g in nd.galaxies:
                    plt.scatter(g.x,g.y,color='red')
            
        plt.show()
    
    
    def plot3D(self):
        """
        3D plot of the XYZ space for the entire tree, showing galaxies
        as red dots. The graph is able to be rotated freely.
        """
        
        from mpl_toolkits.mplot3d import Axes3D
        
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        
        ax.set_xlim(-1,11)
        ax.set_ylim(-1,11)
        ax.set_zlim(-1,11)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        nodeList = self.toList()
        
        
        for nd in nodeList:
            
            if nd.gals > 0:
                for g in nd.galaxies:
                    ax.plot([g.x],[g.y],[g.z],color='red',marker='o',markersize='1')
    
    
    
    def subDiv(self,parent):
        """
        Method for subdividing a node and assigning new child
        node objects to the parent node based on its dimensions
        
        Args:
            parent: expects a pointer to a parent node object
        """
        
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
    
    
    
    
    def evolve(self,step):
        """
        Method for evolving a tree across a time step
        
        Args:
            step: time step to be used for the evolution
            
        Returns:
            returns a new tree that has been populated with evolved galaxies
        """
        
        #Creates a new tree
        newTree = tree(self.xmin,self.xmax,self.ymin,self.ymax,self.zmin,self.zmax)
        
        #Converts time step from per second to per year
        step = step * 31536000
        
        for g in self.galaxies:
            #Calculate acceleration
            g.calcAcceleration(self.root)
            
            #Update position
            newx = 2*g.x - g.oldx + step**2*g.xacc
            newy = 2*g.y - g.oldy + step**2*g.yacc
            newz = 2*g.z - g.oldz + step**2*g.zacc
            
            newTree.addGal(galaxy(newx,newy,newz,g.x,g.y,g.z,g.mass),newTree.root)
        
        newTree.root.calc_COM()
        return newTree
        
        



class galaxy:
    
    def __init__(self,x,y,z,oldx,oldy,oldz,mass):
        """
        Initialization method for creating a galaxy object
        
        Args:
            x: x position of the galaxy
            y: y position of the galaxy
            z: z position of the galaxy
            oldx: x potision of the galaxy in the previous frame
            oldy: y position of the galaxy i nthe previous frame
            oldz: z position of the galaxy in the previous frame
            mass: mass of the galaxy, assumed constant for all galaxies
        """
        
        self.x = x
        self.y = y
        self.z = z
        self.oldx = oldx
        self.oldy = oldy
        self.oldz = oldz
        self.mass = mass
        self.accel = 0
        self.xacc=[]
        self.yacc=[]
        self.zacc=[]
        
    
    def calcAcceleration(self,node,lim=1):
        """
        Method for calculating the acceleration on a galaxy using the Barnes-Hut method
        
        Args:
            node: expects a pointer to a node object, initially the tree root
            lim: default 1, limit for the critical angle below which a node's 
            COM is considered instead of checking the node's children
        """
        
        G = 4.5172e-48
        e = 0.005
        
        #Not a leaf
        if node.children > 0:
            for ch in node.child:
                
                #Check for extraneous cases
                if ch.xcom == None or ch.xcom == 0 or self in ch.galaxies:
                    continue
                
                r = dist(self.x,self.y,self.z,ch.xcom,ch.ycom,ch.zcom)
                
                #Critical theta check
                if r/(ch.xmax-ch.xmin) > lim:
                    self.calcAcceleration(ch)
                else:
                    if ch.contains(self):
                        #Force softening if the node contains the galaxy
                        self.xacc.append(G*ch.mtot*(ch.xcom - self.x) / (r*(r**2+e**2)))
                        self.yacc.append(G*ch.mtot*(ch.ycom - self.y) / (r*(r**2+e**2)))
                        self.zacc.append(G*ch.mtot*(ch.zcom - self.z) / (r*(r**2+e**2)))
                    else:
                        #Regular calculation if the node does not contain the galaxy
                        self.xacc.append(G*ch.mtot*(ch.xcom - self.x) / r**3)
                        self.yacc.append(G*ch.mtot*(ch.ycom - self.y) / r**3)
                        self.zacc.append(G*ch.mtot*(ch.zcom - self.z) / r**3)
        
        #Leaf
        else:
            if self in node.galaxies:
                #Don't calculate anything if the only galaxy in the node is the galaxy being checked
                return
            
            r = dist(self.x,self.y,self.z,node.xcom,node.ycom,node.zcom)
            
            if r < 0.01:
                #Force softened calculation for close galaxies
                self.xacc.append(G*node.mtot*(node.xcom - self.x) / (r*(r**2+e**2)))
                self.yacc.append(G*node.mtot*(node.ycom - self.y) / (r*(r**2+e**2)))
                self.zacc.append(G*node.mtot*(node.zcom - self.z) / (r*(r**2+e**2)))
            
            else:
                #End of the line, use the single galaxies instead of the node
                self.xacc.append(G*node.mtot*(node.xcom - self.x) / r**3)
                self.yacc.append(G*node.mtot*(node.ycom - self.y) / r**3)
                self.zacc.append(G*node.mtot*(node.zcom - self.z) / r**3)
        
        if node.parent == 'Me':
            self.xacc = sum(self.xacc)
            self.yacc = sum(self.yacc)
            self.zacc = sum(self.zacc)


def dist(x1,y1,z1,x2,y2,z2):
    #Simple distance formula method
    return np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)



def save():
    #Imports
    import pickle
    
    #Creates a pickle file with all of the saved instances
    with open(f"series.pk1", 'wb') as output:
       pickle.dump(series, output, pickle.HIGHEST_PROTOCOL)
    with open(f"time.pk1", 'wb') as output:
       pickle.dump(time, output, pickle.HIGHEST_PROTOCOL)
       

def load():
    #Imports
    import pickle
    global series
    global time
    
    with open(f"series",'rb') as input:
        series = pickle.load(input)
    with open(f"time",'rb') as input:
        time = pickle.load(input)


def trackPlot(index=[0]):
    """
    Method for creating a plot that tracks individual galaxies across the time series
    """
    
    from mpl_toolkits.mplot3d import Axes3D
        
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    
    ax.set_xlim(-1,11)
    ax.set_ylim(-1,11)
    ax.set_zlim(-1,11)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    for i in index:
        x = []
        y = []
        z = []
        for s in series:
            x.append(s.galaxies[i].x)
            y.append(s.galaxies[i].y)
            z.append(s.galaxies[i].z)
        
        ax.plot(x,y,z,marker='o')
        





#Create clusters
series = []
cl1 = np.load('galaxies0.npy')
cl2 = np.load('galaxies1.npy')

#Box size
xmin=np.floor(sorted(cl1[:,0])[0])
xmax=np.ceil(sorted(cl1[:,0])[-1])
ymin=np.floor(sorted(cl1[:,1])[0])
ymax=np.ceil(sorted(cl1[:,1])[-1])
zmin=np.floor(sorted(cl1[:,2])[0])
zmax=np.ceil(sorted(cl1[:,2])[-1])

#Constant mass
mass = 10**12

#Timestep
dt = 1e7
time = np.arange(0,10e9,dt)

ref = tree(xmin,xmax,ymin,ymax,zmin,zmax)
series = [tree(xmin,xmax,ymin,ymax,zmin,zmax)]

#Initialize the time series from the data files
for i in range(len(cl1)):
    series[0].addGal(galaxy(cl2[i][0],cl2[i][1],cl2[i][2],cl1[i][0],cl1[i][1],cl1[i][2],mass),series[0].root)
series[0].root.calc_COM()

#Diagnostic checks from initializing the tree
print(f"Tree node depth from 0(root): {series[0].maxDepth}")
print(f"Number of galaxies: {series[0].root.galCount()}")
cList = series[0].toList()
print(f"Number of nodes: {len(cList)}")
print(f"Total COM: {series[0].root.calc_COM()}")

#Initial 3D plot of the cluster
series[0].plot3D()
print(len(series))

#Line-clearing print messages for live updates to code status
from sys import stdout
from time import sleep

print("==========Evolving==========")
for t in time:
    series.append(series[-1].evolve(dt))
    stdout.write(f"\rCurrent Timestep: {t}")
    stdout.flush()
    sleep(0.1)


trackPlot(index=[0,100,-1])