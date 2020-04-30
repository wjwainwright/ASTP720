# -*- coding: utf-8 -*-


class galaxy():
    def __init__(self,name,radius):
        self.name = name
        self.radius = radius
        self.mass = None
        self.volume = None
    
    def calcDensity(self):
        try:
            assert not self.mass == None
        except:
            print(f"Mass has not been calculated for {self.name}")
            return
        try:
            assert not self.volume == None
        except:
            print(f"Volume has not been calculated for {self.name}")
            return
        
        self.density = self.mass/self.volume
        return self.density
    
class spiralGalaxy(galaxy):
    def __init__(self,name,radius,thickness,vrot,Nspiral):
        super().__init__(name,radius)
        self.thickness = thickness
        self.vrot = vrot
        self.Nspiral = Nspiral
    
    def getNspiral(self):
        return self.Nspiral
    
    def calcMass(self):
        G = 6.67e-11
        
        self.mass = self.vrot**2*self.radius/G
        
        return self.mass
    
    def calcVolume(self):
        import numpy as np
        
        self.volume = np.pi*self.radius**2*self.thickness
        
        return self.volume

class ellipticalGalaxy(galaxy):
    def __init__(self,name,radius,velocityDispersion,eccentricity,k=3):
        #Note: radius for elliptical is assumed to be semi-major axis length a
        super().__init__(name,radius)
        self.velocityDispersion = velocityDispersion
        self.k = k
        self.eccentricity = eccentricity
    
    def calcMass(self):
        G = 6.67e-11
        
        self.mass = self.k*self.radius*self.velocityDispersion**2/G
        
        return self.mass
    
    
    def calcVolume(self):
        import numpy as np
        
        a = self.radius
        b = self.radius*np.sqrt(1-self.eccentricity**2)
        
        #Compromise between two types of ellipsoid volumes, a^2b and ab^2
        c = (a+b)/2
        
        self.volume = 4/3*np.pi*a*b*c
        
        return self.volume
 

class dwarfSpheroidal(ellipticalGalaxy):
    pass



#Simulated Galaxy Data to create a density profile
import numpy as np
import matplotlib.pyplot as plt

Ns = 150
Ne = 200
Nds = 15

spirals = []
ellipticals = []
dwarves = []

#Generate Spirals
for i in range(Ns):
    g = spiralGalaxy(f'Spiral_{i:03d}',np.random.uniform(1e19,10e20),np.random.uniform(1e18,10e18),np.random.uniform(400000,700000),np.random.randint(3,10))
    g.calcMass()
    g.calcVolume()
    g.calcDensity()
    spirals.append([g,np.random.uniform(-2e20,2e20),np.random.uniform(-2e20,2e20)])
    assert not g.density == None

#Generate Ellipticals
for i in range(Ne):
    g = ellipticalGalaxy(f'Elliptical_{i:03d}',np.random.uniform(1e17,10e18),np.random.uniform(250000,500000),np.random.uniform(0,0.9))
    g.calcMass()
    g.calcVolume()
    g.calcDensity()
    ellipticals.append([g,np.random.uniform(-2e20,2e20),np.random.uniform(-2e20,2e20)])
    assert not g.density == None
    
#Generate Dwarf Spheroidals
for i in range(Nds):
    g = dwarfSpheroidal(f'Dwarf_{i:03d}',np.random.uniform(1e15,1e17),np.random.uniform(50000,250000),np.random.uniform(0,0.9))
    g.calcMass()
    g.calcVolume()
    g.calcDensity()
    dwarves.append([g,np.random.uniform(-2e20,2e20),np.random.uniform(-2e20,2e20)])
    assert not g.density == None


plt.figure()
plt.title('Positions of the galaxies')
plt.scatter([a[1] for a in spirals],[a[2] for a in spirals],label='Spirals')
plt.scatter([a[1] for a in ellipticals],[a[2] for a in ellipticals],label='Ellipticals',color='red')
plt.scatter([a[1] for a in dwarves],[a[2] for a in dwarves],label='Dwarf Spheroidals',color='orange')
plt.legend()


plt.figure()
plt.title('Positions of the galaxies')
plt.scatter([a[1] for a in spirals],[a[2] for a in spirals],label='Spirals',c=[a[0].density for a in spirals])
plt.scatter([a[1] for a in ellipticals],[a[2] for a in ellipticals],label='Ellipticals',c=[a[0].density for a in ellipticals])
plt.scatter([a[1] for a in dwarves],[a[2] for a in dwarves],label='Dwarf Spheroidals',c=[a[0].density for a in dwarves])
plt.set_cmap('brg')
clb = plt.colorbar()
clb.ax.set_title("Density")









