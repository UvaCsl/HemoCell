# -*- coding: utf-8 -*-
"""
@author: Gabor Zavodszky

Creating a hexagonal patch of material piece.

<-     1 _____ 6     ->
        /     \
<-   2 /  0.   \ 5   ->
       \       /
<-    3 \_____/ 4    ->

"""

import numpy as np

# Mesh and integration defaults
dx = 1e-6
rho = 1025
m = rho * dx**3 # To get similar characteristics to IBM in LBM

class Node (object):
    _id = 0
    def __init__(self, pos=[0.0,0.0,0.0]):
        self.position = np.array(pos)
        self.force = np.array([0.0,0.0,0.0])
        self.velocity = np.array([0.0,0.0,0.0])
        self._fixed = False
        self.Id = Node._id
        Node._id = Node._id+1
        
    def advance(self, dt, mass = m):
        
        if(self._fixed):
            return
            
        #  Euler integration
        self.velocity += self.force * dt / m
        self.position += self.velocity * dt
    
    def resetForce(self):
        self.force.fill(0.0)
        
class Mesh(object):
    
    def __init__(self, l=1.0):
        self.nodes = []
        self.edges = []
        self.faces = []
        self.l_eq = l
        self.createHexPatch()
        
        
    def createHexPatch(self):
        self.nodes = []
        self.edges = []
        self.faces = []
        l_eq = self.l_eq
        dl = l_eq*np.cos(np.pi/6.0)
        
        # Construct nodes
        self.nodes.append(Node([0.0, 0.0, 0.0])) # Origin    
        self.nodes.append(Node([-l_eq/2.0, dl, 0.0]))
        self.nodes.append(Node([-l_eq, 0.0, 0.0]))
        self.nodes.append(Node([-l_eq/2.0, -dl, 0.0]))
        self.nodes.append(Node([l_eq/2.0, -dl, 0.0]))
        self.nodes.append(Node([l_eq, 0.0, 0.0]))
        self.nodes.append(Node([l_eq/2.0, dl, 0.0]))
        
        # Construct edges
        self.edges.append([self.nodes[0],self.nodes[1]])
        self.edges.append([self.nodes[0],self.nodes[2]])
        self.edges.append([self.nodes[1],self.nodes[2]])
        
        self.edges.append([self.nodes[0],self.nodes[3]])
        self.edges.append([self.nodes[2],self.nodes[3]])
        
        self.edges.append([self.nodes[0],self.nodes[4]])
        self.edges.append([self.nodes[3],self.nodes[4]])
        
        self.edges.append([self.nodes[0],self.nodes[5]])
        self.edges.append([self.nodes[4],self.nodes[5]])
        
        self.edges.append([self.nodes[0],self.nodes[6]])
        self.edges.append([self.nodes[5],self.nodes[6]])
        
        self.edges.append([self.nodes[6],self.nodes[1]])


        # Construct faces
        self.faces.append([self.nodes[0],self.nodes[1],self.nodes[2]])
        self.faces.append([self.nodes[0],self.nodes[2],self.nodes[3]])
        self.faces.append([self.nodes[0],self.nodes[3],self.nodes[4]])
        self.faces.append([self.nodes[0],self.nodes[4],self.nodes[5]])
        self.faces.append([self.nodes[0],self.nodes[5],self.nodes[6]])
        self.faces.append([self.nodes[0],self.nodes[6],self.nodes[1]])
        

    def getVerticalStretch(self):
        d = self.nodes[5].position - self.nodes[2].position
        return np.sqrt(np.inner(d,d)) 
    
    def calcFaceArea(self, v0, v1, v2):
        cr = np.cross(v1-v0, v2-v0) 
        return 0.5 * np.sqrt(np.inner(cr,cr))
    
    def calcTotalArea(self):
        a = 0.0
        for f in self.faces:
            a += self.calcFaceArea(f[0], f[1], f[2])
        return a
    
    def calcFaceNormal(self, v0, v1, v2):
        cr = np.cross(v1-v0, v2-v0)
        return cr / np.sqrt(np.inner(cr,cr))
        
    def advance(self, dt):
        for n in self.nodes:
            n.advance(dt)
            
    def resetForces(self):
        for n in self.nodes:
            n.resetForce()
            
    def getAreaEq(self):
        # Assuming equilateral triangles
        return np.sqrt(3.0)/4.0 * self.l_eq * self.l_eq
    
    def applyStretchForce(self, force):
        force_pp = force / 3.0
        self.nodes[1].force += [-force_pp, 0, 0]
        self.nodes[2].force += [-force_pp, 0, 0]
        self.nodes[3].force += [-force_pp, 0, 0]
        self.nodes[4].force += [ force_pp, 0, 0]
        self.nodes[5].force += [ force_pp, 0, 0]
        self.nodes[6].force += [ force_pp, 0, 0]

    def plot(self):
        x, y = self.getPositions()
        import matplotlib.pyplot as plt
        plt.plot(x,y, 'o')
        
    def getPositions(self):
        x=[]; y=[]
        for i in range(len(self.nodes)):
            x.append(self.nodes[i].position[0])
            y.append(self.nodes[i].position[1])
        return (x,y)