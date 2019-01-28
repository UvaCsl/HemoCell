import numpy as np
from numpy.linalg import eig, inv
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import math
import pdb


def get_mario_curves(what="full",datapath = "/Users/bczaja/Work/PHD/michigan/EL_curves/"):
    if what=="full":
        data = np.loadtxt(datapath +"Ekcta_100.csv","float", delimiter=";", skiprows=1)
    
    return data

from numpy.linalg import eig, inv

def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a

def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])


def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))


def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])

def ellipse_angle_of_rotation2( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    if b == 0:
        if a > c:
            return 0
        else:
            return np.pi/2
    else: 
        if a > c:
            return np.arctan(2*b/(a-c))/2
        else:
            return np.pi/2 + np.arctan(2*b/(a-c))/2
        
def ellipsePlots(center, phi, axes):
    u=center[0]       #x-position of the center
    v=center[1]      #y-position of the center
    a=axes[0]      #radius on the x-axis
    b=axes[1]    #radius on the y-axis
    t_rot=phi #rotation angle

    t = np.linspace(0, 2*np.pi, 1000)
    Ell = np.array([a*np.cos(t) , b*np.sin(t)])  
         #u,v removed to keep the same center location
    R_rot = np.array([[np.cos(t_rot) , -np.sin(t_rot)],[np.sin(t_rot) , np.cos(t_rot)]])  
         #2-D rotation matrix

    Ell_rot = np.zeros((2,Ell.shape[1]))
    for i in range(Ell.shape[1]):
        Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])

    return u+Ell_rot[0,:], v+Ell_rot[1,:] 


def elongation_index(x,y,dx=0.5e-6):
    new =[]
    slice_dx =dx
    for i in np.arange(np.min(x),np.max(x)+slice_dx,slice_dx):
        #Get the indexes of the masked region
        idx_mask = np.where((x <= i + slice_dx) & (x> i))
        #Make sure that data exists in this range
        if len(y[idx_mask]) >0:
        
        
            
            #Get the unique min and max
            ymax = np.unique(np.max(y[idx_mask]))[0]
            ymin = np.unique(np.min(y[idx_mask]))[0]
                    
            #Find again the index of the min and max
            final_idx_max = np.where((y == ymax) & (x<=i+slice_dx) & (x>i))[0][0]
            final_idx_min = np.where((y == ymin) & (x<=i+slice_dx) & (x>i))[0][0]
            
            #get the corresponding xcoordinate for the ymax
            x_for_ymax = x[final_idx_max]
            x_for_ymin = x[final_idx_min]
    
            new.append(np.array([x_for_ymax,ymax]))
            new.append(np.array([x_for_ymin,ymin]))
        
    new = np.array(new)
    
    
    a = fitEllipse(new[:,0],new[:,1])
    center = ellipse_center(a)
    phi = ellipse_angle_of_rotation2(a)
    axes = ellipse_axis_length(a)
    xF,yF = ellipsePlots(center,phi,axes)
    
    major = axes[0]
    minor = axes[1]
    El = (major - minor)/(major+minor)
    
    return(major,minor,El)





def elongation_index_OLD(x,y,z,plot=False):
    '''returns elongation index, major axis, and minor axis'''
    
    #if plot:
    #    plt.close()
    #    plt.scatter(x,y)
    #    plt.savefig("/Users/bczaja/Desktop/elongation_index_snapshot_crap.png",dpi=200)
    #    plt.show()

    #Center the coordenates at the center of the cell
    x -= np.mean(x)
    y -= np.mean(y)
    z -= np.mean(z)
    
    #Find the farthers two points via scipy
    XY = np.column_stack((x,y))
    D=pdist(XY)
    D=squareform(D)
    temp = np.where(D == D.max())
    #x & y coordinates of the farthers points
    xfar = x[temp[0]]
    yfar = y[temp[0]]
    yfit = np.poly1d(np.polyfit(xfar,yfar,1))
    xp = np.linspace(-100,100,1000)
    
    idx_origin = yfit(xp) == np.min(np.abs(yfit(xp)))
    if not yfit(xp)[idx_origin]:
        idx_origin = yfit(xp) == -np.min(np.abs(yfit(xp)))
    
    hypotenuse = np.sqrt(np.max(xfar)**2 + np.max(yfar)**2)
    adjacent = np.max(xfar)
    theta =   math.acos(adjacent/hypotenuse)
    
    #Rotate the coordinates
    x_prime = x*math.cos(theta) + y*math.sin(theta)
    y_prime = -1.*x*math.sin(theta) + y*math.cos(theta)
    #Also the farthest points    
    xfar_prime = xfar*math.cos(theta) + yfar*math.sin(theta)
    yfar_prime = -1.*xfar*math.sin(theta) + yfar*math.cos(theta)
    
    #now find the major axis
    A = np.sqrt(np.max(xfar_prime)**2 + np.max(yfar_prime)**2) + np.sqrt(np.min(xfar_prime)**2 + np.min(yfar_prime)**2)
    
    #Find the minor Axis
    #Scan the length for the longest axis
    B = 0
    for i in x_prime:
        slice_mask = (x_prime <= i + 1) & (x_prime >= i)
        B_top =    np.max(y_prime[slice_mask])
        B_bottom = np.min(y_prime[slice_mask])
        
        if (np.sqrt(B_top**2 + B_bottom**2) > B):
            B = np.sqrt(B_top**2 + B_bottom**2)

    
    #Calculate Elongation Index
    elongation_index = (A-B)/(A+B)
    
    
    #Bad method1
    #Find the midpoint along the x axis
    #xmid = np.mean(np.array([np.max(xfar_prime),np.min(xfar_prime)]))
    #xmid_mask = (x_prime <= xmid +1) & (x_prime >= xmid -1)
    #stripex = x_prime[xmid_mask]
    #stripey = y_prime[xmid_mask]
    #B_top = np.max(y_prime[xmid_mask])
    #B_bottom = np.min(y_prime[xmid_mask])
    #B2 = np.sqrt(B_top**2 + B_bottom**2)
    
    #Bad method2
    x_minor_top = x_prime[y_prime == np.max(y_prime)][0]
    y_minor_top = y_prime[y_prime == np.max(y_prime)][0]
    x_minor_bottom = x_prime[y_prime == np.min(y_prime)][0]
    y_minor_bottom = y_prime[y_prime == np.min(y_prime)][0]
    #B = (np.sqrt(x_minor_top**2 + y_minor_top**2) + np.sqrt(x_minor_bottom**2 + y_minor_bottom**2))
        

    
    
    if plot:
        plt.close()
        plt.scatter(x_prime,y_prime)
        plt.scatter(x_minor_top,y_minor_top,c='purple')
        plt.scatter(x_minor_bottom,y_minor_bottom,c='purple',label ="B = %.6f"%(B))
        plt.scatter(stripex,stripey,c='cyan',label ="B = %.6f"%(B))
        plt.scatter(xfar_prime,yfar_prime,c = 'r',label ="A = %.6f"%(A))
        plt.xlim(-22,22)
        plt.ylim(-22,22)
        plt.legend()
        plt.title("Elongation index = %.3f"%(elongation_index))
        plt.savefig("/Users/bczaja/Desktop/elongation_index_snapshot.png",dpi=200)
        plt.show()

    return (A,B,elongation_index)










