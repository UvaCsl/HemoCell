# -*- coding: utf-8 -*-
"""
@author: Gabor Zavodszky
"""
import numpy as np
from getModuli.cellModel import CellModel
from getModuli.mesh import Mesh

class Cell(CellModel):
    
    def __init__(self, l=1.0e-6):
        self.mesh = Mesh(l=l)
        self.l_eq = self.mesh.l_eq
        self.area_eq = self.mesh.getAreaEq()
        self.k_link = 15.0 * self.kBT / self.persistenceLengthFine
        self.k_area = 5.0 * self.kBT / self.l_eq**2
        self.eta_m = 1e-9   

        
    def calcConstitutiveForces(self):
        #nodes = self.mesh.nodes
        edges = self.mesh.edges
        faces = self.mesh.faces
           
        # Edge forces
        for e in range(len(edges)):
            v0 = edges[e][0].position
            v1 = edges[e][1].position
            edge_vec = v1 - v0
            edge_len = np.sqrt(np.inner(edge_vec,edge_vec)) #norm
            edge_dir = edge_vec/edge_len
            
            edge_frac = (edge_len - self.l_eq)/self.l_eq
            edge_force = edge_dir * self.k_link * (edge_frac + edge_frac/(9.0-edge_frac**2))
            edges[e][0].force += edge_force
            edges[e][1].force -= edge_force
                 
            # Membrane viscosity
            v_rel = edges[e][1].velocity - edges[e][0].velocity
            v_edge_rel = np.dot(v_rel, edge_dir)
            edge_visc_f = self.eta_m * v_edge_rel * edge_dir
            edges[e][0].force += edge_visc_f
            edges[e][1].force -= edge_visc_f    
        
        # Area forces
        for f in range(len(faces)):
            v0 = faces[f][0].position
            v1 = faces[f][1].position
            v2 = faces[f][2].position
            
            face_area = self.mesh.calcFaceArea(v0, v1, v2)
            area_frac = (face_area-self.area_eq)/self.area_eq
            force_mag = self.k_area * (area_frac + area_frac/(0.09-area_frac**2))
            centroid = [ (v0[0]+v1[0]+v2[0])/3.0,
                         (v0[1]+v1[1]+v2[1])/3.0,
                         (v0[2]+v1[2]+v2[2])/3.0]
            av0 = centroid - v0
            av1 = centroid - v1
            av2 = centroid - v2
            
            faces[f][0].force += force_mag * av0
            faces[f][1].force += force_mag * av1
            faces[f][2].force += force_mag * av2
    
    def advance(self, dt):
        self.mesh.advance(dt)