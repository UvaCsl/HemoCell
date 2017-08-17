# -*- coding: utf-8 -*-
"""
@author: Gabor Zavodszky
"""
import numpy as np
from getModuli.rbcHO import RbcHO
# from getModuli.wbcHO import WbcHO 

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
    # model = WbcHo(l_eq)
    
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
    
    stretch_end = model.mesh.getVerticalStretch()
    dfl = 2.0 * l_eq#2.0 * l_eq*np.cos(np.pi/6.0)
    dl = stretch_end - 2.0 * l_eq                                           
    young = force * 1e6 * 2.0*l_eq / (dfl * dl)
        
    print "After ", iteration," iterations -> final stretch: ", stretch_end, " (", stretch_end / (2.0*l_eq) * 100.0 ,"%), force: ",  f_stretch    
    print "Approx. Young mod. [uN/m]: ", young   
    return (timeStep, stretch, model.mesh.getPositions())


def static_stretch(disloc):
    
    # Constructing patch
    model = RbcHO(l_eq)
    # model = WbcHo(l_eq)
    a0 = model.mesh.calcTotalArea()
    
    # Move verticies by disloc
#    
#    # Left side
#    for i in range(1,4):
#        model.mesh.nodes[i].position[0] -= disloc
#    # Right side
#    for i in range(4,7):
#        model.mesh.nodes[i].position[0] += disloc

    dl = disloc * np.cos(np.pi/6.0)
    dd = disloc * 0.5
    model.mesh.nodes[1].position += [-dd,  dl, 0.0]
    model.mesh.nodes[2].position += [-disloc, 0.0, 0.0]
    model.mesh.nodes[3].position += [-dd, -dl, 0.0]
    model.mesh.nodes[4].position += [dd, -dl, 0.0]
    model.mesh.nodes[5].position += [disloc, 0.0, 0.0]
    model.mesh.nodes[6].position += [dd, dl, 0.0]   

    a = model.mesh.calcTotalArea()
                                          
    # Compute arising forces
    model.mesh.resetForces()
    model.calcConstitutiveForces()
    
    force = [0.0, 0.0, 0.0]
    for i in range(1,4):
        force += model.mesh.nodes[i].force

    force_mag = np.sqrt(np.inner(force,force))
    
    dfl = 2.0 * l_eq*np.cos(np.pi/6.0) #2.0 * l_eq
    da = a - a0
    compr2d = force_mag * 1e6 * a0 / (dfl * da)
    
    print "Approx. 2D compression mod. [uN/m]: ", compr2d 
    
#    stretch_end = model.mesh.getVerticalStretch()                             
#    dl = stretch_end - 2.0 * l_eq                                           
#    compr = force_mag * 1e6 * 2.0*l_eq / (2.0 * l_eq * dl)  
#    print "Force resonse: ", force_mag, ", for stretch: ",  stretch_end
#    print "Approx. lin. compression mod. [uN/m]: ", compr
                                   
                              
    return (compr2d, force, model.mesh.getPositions())      


def static_shear(disloc):
    
    # Constructing patch
    model = RbcHO(l_eq)
    # model = WbcHo(l_eq)
    
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
    shear = force_mag * 1e6 * dy / (2.0 * l_eq * disloc)
    
    print "Force resonse: ", force_mag, ", for shear: ",  dy/disloc
    print "Approx. shear mod. [uN/m]: ", shear
                              
    return (shear, force, model.mesh.getPositions())                      

if __name__ == "__main__":
    
    print "\n ### Test: Dynamic stretching ###"
    f_stretch = 3.0e-12  # [N]
    t, s, p = dynamic_stretch(f_stretch)

    
    print "\n ### Test: Static stretching ###"
    K, f, p2 = static_stretch(l_eq*0.1)
    
    print "\n ### Test: Static shearing ###"
    mu0, f3, p3 = static_shear(l_eq*0.1)
    
    young = 4.0 * K * mu0 / (K + mu0)
    poisson = (3.0 * K - 2.0 * mu0) / (2.0*(3.0*K + mu0))
    
    print "\n ### Derived quantities ###"
    print "Young modulus: ", young, " [uN/m]"
    print "Poisson ratio: ", poisson
    
    # PLotting
    # import matplotlib.pyplot as plt    
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # ax1.scatter(p[0], p[1], c='b', marker="o", label='dyn. stretch')
    # ax1.scatter(p2[0], p2[1], c='r', marker="o", label='stat. stretch')
    # ax1.scatter(p3[0], p3[1], c='g', marker="o", label='stat. shear')

    # # Ref.
    # model = RbcHO(l_eq)
    # # model = WbcHo(l_eq)
    # pr = model.mesh.getPositions()
    # ax1.scatter(pr[0], pr[1], c='y', marker="v", label='Undeformed')

    # plt.xlim([-2.0*l_eq, 2.0*l_eq])
    # plt.ylim([-2.0*l_eq, 2.0*l_eq])
    # plt.legend(loc='upper left');
    # plt.axes().set_aspect('equal', 'datalim')
    # plt.show()
