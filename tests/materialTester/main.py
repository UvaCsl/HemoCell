# -*- coding: utf-8 -*-
"""
@author: Gabor Zavodszky
"""
import numpy as np
from getModuli.rbcHO import RbcHO 

l_eq = 0.5e-6


def dynamic_stretch(force):
  
    # Parameters
    v_eps = 1e-10   # stopping creterion
    dt = 5e-8       # integration time-step  
    
    stretch = []
    timeStep = []
    iteration = 0
    
    # Constructing patch
    model = RbcHO(l_eq)
    
    # Starting simulation
    while True:
    #for i in range(4000):

        if iteration % 2000 == 0:
            timeStep.append(iteration*dt)
            stretch.append(model.mesh.getVerticalStretch() * 0.5)
            print "Time: ", timeStep[-1], "s, stretch: ", stretch[-1]
        
        iteration += 1
        
        model.mesh.resetForces()
        model.mesh.applyStretchForce(force)
        model.calcConstitutiveForces()
        model.advance(dt)
        
        # Check stop criterion
        v = model.mesh.nodes[2].velocity
        if np.sqrt(np.inner(v,v))<v_eps:
            break
    
    stretch_end = model.mesh.getVerticalStretch() * 0.5
    
    young = f_stretch * 1e6 * l_eq / (l_eq * (stretch_end - l_eq))
        
    print "After ", iteration," iterations -> final stretch: ", stretch_end, " (", stretch_end / l_eq * 100.0 ,"%), force: ",  f_stretch    
    print "Approx. Young mod. [uN/m]: ", young   
    return (timeStep, stretch, model.mesh.getPositions())


def static_stretch(disloc):
    
    # Constructing patch
    model = RbcHO(l_eq)
    
    # Move verticies by disloc
    
    # Left side
    for i in range(1,4):
        model.mesh.nodes[i].position[0] -= disloc
    # Right side
    for i in range(4,7):
        model.mesh.nodes[i].position[0] += disloc
                                          
    # Compute arising forces
    model.mesh.resetForces()
    model.calcConstitutiveForces()
    
    force = [0.0, 0.0, 0.0]
    for i in range(1,4):
        force += model.mesh.nodes[i].force

    stretch_end = model.mesh.getVerticalStretch() * 0.5    
                             
    force_mag = np.sqrt(np.inner(force,force))
    compr = force_mag * 1e6 * l_eq / (l_eq * (stretch_end - l_eq))
    
    print "Force resonse: ", force_mag, ", for stretch: ",  stretch_end
    print "Approx. compression mod. [uN/m]: ", compr
                              
    return (young, force, model.mesh.getPositions())      


def static_shear(disloc):
    
    # Constructing patch
    model = RbcHO(l_eq)
    
    # Move verticies a little
    
    # Middle row
    for i in [2,0,5]:
        model.mesh.nodes[i].position[0] += disloc * 0.5
    # Top row
    for i in [1,6]:
        model.mesh.nodes[i].position[0] += disloc

    dy = 2.0 * l_eq*np.cos(np.pi/6.0)    
                                      
    # Compute arising forces
    model.mesh.resetForces()
    model.calcConstitutiveForces()
    
    force = [0.0, 0.0, 0.0]
    for i in [1,6,5]:  
        force += model.mesh.nodes[i].force
                                 
    force_mag = np.sqrt(np.inner(force,force))
    shear = force_mag * 1e6 * dy / (dy * disloc)
    
    print "Force resonse: ", force_mag, ", for shear: ",  dy/disloc
    print "Approx. shear mod. [uN/m]: ", shear
                              
    return (shear, force, model.mesh.getPositions())                      

if __name__ == "__main__":
    
    print "\n ### Test: Dynamic stretching ###"
    f_stretch = 3.0e-12  # [N]
    t, s, p = dynamic_stretch(f_stretch)

    
    print "\n ### Test: Static stretching ###"
    young, f, p2 = static_stretch(l_eq*0.1)
    
    print "\n ### Test: Static shearing ###"
    shear, f3, p3 = static_shear(l_eq*0.1)
    
    # PLotting
    import matplotlib.pyplot as plt    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(p[0], p[1], c='b', marker="o", label='dyn. stretch')
    ax1.scatter(p2[0], p2[1], c='r', marker="o", label='stat. stretch')
    ax1.scatter(p3[0], p3[1], c='g', marker="o", label='stat. shear')

    # Ref.
    model = RbcHO(l_eq)
    pr = model.mesh.getPositions()
    ax1.scatter(pr[0], pr[1], c='y', marker="v", label='Undeformed')

    plt.xlim([-2.0*l_eq, 2.0*l_eq])
    plt.ylim([-2.0*l_eq, 2.0*l_eq])
    plt.legend(loc='upper left');
    plt.show()
