
'''
============================================================================
Helper functions for the analysis file Run.py
Authors: Benjamin Gorissen
Harvard University
email: bgorissen@seas.harvard.edu
Last update: Sree Patiballa

Tested on Abaqus 2019
============================================================================
'''

## Function for reading CSV file
def csv2list(filename):
  with open(filename,'r') as csvfile:
    csvreader=csv.reader(csvfile)
    data = list(csvreader)
  return data

def plane_to_cyl(x,radius):
    theta = x/radius    
    x_new = radius*cos(theta)
    zDisp = radius*sin(theta)
    xDisp = x_new-x
    yDisp = 0.
    return (xDisp, yDisp, zDisp)

def squarepatch(width, height, r,theta,y, patch_angle):
    x_org = r*theta
    y_org = y
    coords_out = []

    x_coords=[width/2.,-width/2.]
    y_coords=[height/2.,-height/2.]
    for x in x_coords:
        for y in y_coords:
            coords_out.append([cos(patch_angle)*x-sin(patch_angle)*y+x_org, sin(patch_angle)*x+cos(patch_angle)*y+y_org])
    coords_out[2], coords_out[3]=coords_out[3], coords_out[2]
    return coords_out

def mean(a):
    return sum(a) / len(a)

def is_patch_inbounds(PATCH,x_max,y_max): #x-coordinates go from 0 to CYL_CIRCUM and y-coordinates go from 0 to CYL_CYL_LENGTH
    output = True
    COORDS = PATCH[0]
    for COORD in COORDS:
        if COORD[0]>x_max or COORD[0]<0.:
            output = False
            break
        if COORD[1]>y_max or COORD[1]<0.:
            output = False
            break
    return output

def divide_patch(PATCH,Radius):
    coords_original = PATCH[0]
    coords_all = []
    for indexx,coords1 in enumerate(coords_original):
        coords2 = coords_original[indexx-1]
        coords_all.append(coords1)
        if (max(coords1[0],coords2[0]) > 2*pi*Radius or min(coords1[0],coords2[0]) < 0):
            newcoords=interpolate_points(coords1,coords2,Radius)
            if newcoords: 
                for newcoord in newcoords:
                    coords_all.append(newcoord)

    newcoords = create_new_patches(coords_all, Radius)
    newpatches = []
    for newc in newcoords:
        newpatches.append([newc, PATCH[-1]])
    return newpatches



def interpolate_points(coords1,coords2,Radius): #Function to create now points along a line between coords1 and coords2 at an interval of 2*pi*Radius
    if coords2[0] < coords1[0]: #sort coords1 and coords2, coords1 should be left of coords 2
        coords_temp = coords1
        coords1 = coords2
        coords2 = coords_temp


    left_x=float(int(coords1[0]/(2*pi*Radius)))*2*pi*Radius
    if coords1[0] > 0:
        left_x = left_x+2*pi*Radius
    right_x = float(int(coords2[0]/(2*pi*Radius)))*2*pi*Radius
    if coords2[0] < 0:
        right_x = right_x-2*pi*Radius

    newcoords = []

    x_coords = np.arange(left_x,right_x+0.001,2.*pi*Radius).tolist()
    y_coords = np.interp(x_coords,[coords1[0],coords2[0]],[coords1[1],coords2[1]])


    for i, xx in enumerate(x_coords):
        newcoords.append([xx,y_coords[i]])
    return(newcoords)

def create_new_patches(coords_all,Radius):
    coords_all = sorted(coords_all, key=lambda x: x[0]) # sort all the values according to their x values
    xmin = coords_all[0][0]
    xmax = coords_all[-1][0]
    xleftborder = round(xmin/(2.*pi*Radius))*2.*pi*Radius
    if xmin<xleftborder:
        xleftborder = xleftborder-2.*pi*Radius
    xrightborder = round(xmax/(2.*pi*Radius))*2.*pi*Radius
    if xmax>xrightborder:
        xrightborder = xrightborder+2.*pi*Radius
    xticks = np.arange(xleftborder,xrightborder+0.01,2.*pi*Radius)
    NEWCOORDS=[]
    for i, xx in enumerate(xticks):
        if i !=0: 
            NEWCOORDSPATCH = []
            for coords in coords_all:
                if coords[0]>(xticks[i-1]-0.001) and coords[0]<(xticks[i]+0.001):
                    coordspatch = copy.copy(coords)
                    NEWCOORDSPATCH.append(coordspatch)
            NEWCOORDS.append(NEWCOORDSPATCH)
    for COORDS in NEWCOORDS:
        COORD_CENTR = map(mean, zip(*COORDS))
        shift = floor(COORD_CENTR[0]/(2*pi*Radius))
        for COORD in COORDS:
            COORD[0] = COORD[0]-shift*2*pi*Radius
    return NEWCOORDS

def make_ccw(PATCH):
    COORDS = PATCH[0]
    COORD_CENTR = map(mean, zip(*COORDS))
    COORDSPOLAR=[]
    for COORD in COORDS:
        xrel = COORD[0]-COORD_CENTR[0]
        yrel = COORD[1]-COORD_CENTR[1]
        COORDPOLAR=[np.sqrt(xrel**2+yrel**2), np.arctan2(yrel,xrel)]
        COORDSPOLAR.append(COORDPOLAR)
    COORDSPOLAR = sorted(COORDSPOLAR, key=lambda x: x[1])
    COORDS = []
    for COORDPOLAR in COORDSPOLAR:
        xrel=COORDPOLAR[0]*np.cos(COORDPOLAR[1])
        yrel=COORDPOLAR[0]*np.sin(COORDPOLAR[1])
        COORD = [xrel+COORD_CENTR[0], yrel+COORD_CENTR[1]]
        COORDS.append(COORD)
    PATCH[0] = COORDS
    return(PATCH)



#-------------------------------------------------------------------------------
# RODRIGUES ROTATION
# - Rotate given points based on a starting and ending vector
# - Axis k and angle of rotation theta given by vectors n0,n1
#   P_rot = P*cos(theta) + (k x P)*sin(theta) + k*<k,P>*(1-cos(theta))
#-------------------------------------------------------------------------------

def rodrigues_rot(P, n0, n1):
    
    # If P is only 1d array (coords of single point), fix it to be matrix
    if P.ndim == 1:
        P = P[np.newaxis,:]
    
    # Get vector of rotation k and angle theta
    n0 = n0/np.linalg.norm(n0)
    n1 = n1/np.linalg.norm(n1)
    k = np.cross(n0,n1)
    k = k/np.linalg.norm(k)
    theta = np.arccos(np.dot(n0,n1))
    
    # Compute rotated points
    P_rot = np.zeros((len(P),3))
    for i in range(len(P)):
        P_rot[i] = P[i]*np.cos(theta) + np.cross(k,P[i])*sin(theta) + k*np.dot(k,P[i])*(1-np.cos(theta))

    return P_rot


#-------------------------------------------------------------------------------
# FIT CIRCLE 2D
# - Find center [xc, yc] and radius r of circle fitting to set of 2D points
# - Optionally specify weights for points
#
# - Implicit circle function:
#   (x-xc)^2 + (y-yc)^2 = r^2
#   (2*xc)*x + (2*yc)*y + (r^2-xc^2-yc^2) = x^2+y^2
#   c[0]*x + c[1]*y + c[2] = x^2+y^2
#
# - Solution by method of least squares:
#   A*c = b, c' = argmin(||A*c - b||^2)
#   A = [x y 1], b = [x^2+y^2]
#-------------------------------------------------------------------------------

def fit_circle_2d(x, y, w=[]):
    
    A = array([x, y, np.ones(len(x))]).T
    b = x**2 + y**2
    
    # Modify A,b for weighted least squares
    if len(w) == len(x):
        W = diag(w)
        A = np.dot(W,A)
        b = np.dot(W,b)
    # Solve by method of least squares
    c = np.linalg.lstsq(A,b)[0]
    
    # Get circle parameters from solution c
    xc = c[0]/2
    yc = c[1]/2
    r = sqrt(c[2] + xc**2 + yc**2)
    return xc, yc, r



#-------------------------------------------------------------------------------------
# Frenet Code
# To compute the torsion and curvature of a centerline coordinates
#-------------------------------------------------------------------------------------
import numpy as np 
from numpy import genfromtxt 

def frenet(x,y,z):

    # FRENET - Frenet-Serret Space Curve Invarients   
    # Code translated from matlab to python by Robert Baines 4/26/2021 

    # INPUT : X Y AND Z  of type np.array 

    # SPEED OF CURVE
    dx = np.gradient(x)
    dy = np.gradient(y)
    dz = np.gradient(z)
    dr = np.vstack((dx, dy, dz))

    # Accelleration OF CURVE
    ddx = np.gradient(dx)
    ddy = np.gradient(dy)
    ddz = np.gradient(dz)
    ddr = np.vstack((ddx, ddy, ddz))

    dddx = np.gradient(ddx)
    dddy = np.gradient(ddy)
    dddz = np.gradient(ddz)
    dddr = np.vstack((dddx, dddy, dddz))

    # TANGENT
    T = np.divide(dr, mag(dr,3))

    # DERIVIATIVE OF TANGENT
    dTx =  np.gradient(T[0,:])
    dTy =  np.gradient(T[1,:])
    dTz =  np.gradient(T[2,:])

    dT = np.vstack((dTx, dTy, dTz))


    # NORMAL
    N = np.divide(dT, mag(dT,3))
    #BINORMAL
    B = np.cross( T,N,axisa=0, axisb=0)

    # CURVATURE
    k =  np.divide(  mag(  np.transpose( np.cross(dr,ddr,axisa=0, axisb=0) ), 1 )   ,   np.power((mag(dr,1)),3)   )  # ROB: validated equation here https://math.libretexts.org/Bookshelves/Calculus/Supplemental_Modules_(Calculus)/Vector_Calculus/2%3A_Vector-Valued_Functions_and_Motion_in_Space/2.3%3A_Curvature_and_Normal_Vectors_of_a_Curve
    
    # TORSION    
    t =   np.divide (   vertdot(  np.transpose(np.cross(dr, ddr,axisa=0, axisb=0))  , dddr)      ,    np.power(  mag(  np.transpose( np.cross(dr,ddr,axisa=0, axisb=0) ), 1 )       , 2 )   ) # ROB: and here https://www.math.upenn.edu/~wziller/math114f13/ch13-5+6.pdf

    return k,t


#################

def mag(T,n):
# MAGNATUDE OF A VECTOR (Nx3)
#  M = mag(U)

    N =    np.power (  np.sum( np.power ( np.absolute(T) , 2 )  ,0)   , 0.5 )

    d = np.argwhere(N==0) # find where we have zero


    N[d] = 2.2204e-16*np.ones(np.shape(d)) # replace zero with a small number so we don't break during division. 

    N = np.tile(N,(n,1)) 

    return N 

#################

def vertdot(A, B):
    #row-wise dot-product of A and B
    N=np.zeros(np.shape(A)[1])
    for i in range ( np.shape(A)[1] ): 

        N[i] = np.dot(A[:,i], B[:,i])
    
    return N 