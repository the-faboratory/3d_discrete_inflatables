
'''
============================================================================
Python code for reading out the results of a simulation of a cylindrical shell with patches
NB: The units are N and mm (MPa)
Authors: Benjamin Gorissen
Harvard University
email: bgorissen@seas.harvard.edu
Last update: Sree Patiballa, Yale University
Current affiliation: University of Alabama, spatiballa@ua.edu

Tested on Abaqus 2019
============================================================================
'''
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from time import sleep
from odbAccess import openOdb

import numpy as np

PROJECT_NAME = 'Cylinder';
JOB_NAME      = PROJECT_NAME+'_job'
STEPNAME_PRESSURE ='Holding-Step'

CYL_RADIUS      = 3.18 # 3. Radius of cylinder in mm
# CYL_LENGTH      = 300.  # Length of cylinder in mm
CYL_LENGTH      = 292.  # Length of cylinder in mm
CYL_THICKNESS   = 0.183  # 0.33 Cylinder thickness in mm
PATCH_THICKNESS = 0.49
CYL_CIRCUM      = 2*pi*CYL_RADIUS
MESHSIZE        = CYL_CIRCUM/8. #Size of the elements
TOLERANCE       = MESHSIZE/4.

NUMOUTNODES   = 64 #Number of central nodes to determine output shape ()

'''
------------------------------------------------------------------------
READ OUT
------------------------------------------------------------------------
'''

execfile('HelperFunctions.py')

odb = session.openOdb('Cylinder_job.odb');

X_CENTERS = []
Y_CENTERS = []
Z_CENTERS = []

def fit_circle_2d(x, y, w=[]):
    A = np.array([x, y, np.ones(len(x))]).T
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


for index in range(NUMOUTNODES):
	NODES = odb.rootAssembly.instances['CYLINDER_ASSEMBLY'].nodeSets['SET-OUTPUT-'+str(index)].nodes
	x=[]
	y=[]
	z=[]
	for NODE in NODES:
		LABEL = NODE.label
		historyRegion = odb.steps[STEPNAME_PRESSURE].historyRegions['Node '+NODE.instanceName+'.'+str(NODE.label)]
		out_his_1 = historyRegion.historyOutputs['U1'].data;
		out_his_2 = historyRegion.historyOutputs['U2'].data;
		out_his_3 = historyRegion.historyOutputs['U3'].data;
		ux = out_his_1[-1][-1]
		uy = out_his_2[-1][-1]
		uz = out_his_3[-1][-1]
		x0,y0,z0 = historyRegion.point.node.coordinates
		xp=(x0+ux)
		yp=(y0+uy)
		zp=(z0+uz)
		x.append(xp)
		y.append(yp)
		z.append(zp)
	#Instead of taking the average of the points, do a circle fit. 
	xarr = np.array(x)
	yarr = np.array(y)
	zarr = np.array(z)
	arraytuple = (xarr,yarr,zarr)
	P=np.transpose(np.vstack(arraytuple))
	P_mean = P.mean(axis=0)
	P_centered = P - P_mean
	U,s,V = np.linalg.svd(P_centered)
	normal = V[2,:]
	d = -np.dot(P_mean, normal)
	P_xy = rodrigues_rot(P_centered, normal, [0,0,1])
	xc, yc, r = fit_circle_2d(P_xy[:,0], P_xy[:,1])
	C = rodrigues_rot(np.array([xc,yc,0]), [0,0,1], normal) + P_mean
	X_CENTERS.append(C[0][0])
	Y_CENTERS.append(C[0][1])
	Z_CENTERS.append(C[0][2])

# Make sure that the y-axis is alligend

xarr = np.array(X_CENTERS)
yarr = np.array(Y_CENTERS)
zarr = np.array(Z_CENTERS)

arraytuple = (xarr,yarr,zarr)
P=np.transpose(np.vstack(arraytuple))

normal_rot = np.array([X_CENTERS[19],Y_CENTERS[19],Z_CENTERS[19]])
print(normal_rot)

# works well for 3d shapes
P_rot = rodrigues_rot(P, normal_rot, [0,1,1])

# works well for 2d  and simple shapes
# P_rot = rodrigues_rot(P, normal_rot, [0,0,1])

print(normal_rot)
print('here')
print(P)
print('rotation')
print(P_rot)

P_rot = np.transpose(P_rot)

X_CENTERS=P_rot[0].tolist()
Y_CENTERS=P_rot[1].tolist()
Z_CENTERS=P_rot[2].tolist()

def_coord = []
for index in range(len(X_CENTERS)):
        def_coord.append([X_CENTERS[index],Y_CENTERS[index],Z_CENTERS[index]])


np.savetxt("centerline_coord.csv", def_coord, delimiter=",")



'''
------------------------------------------------------------------------
Write resutls
------------------------------------------------------------------------
'''
# f=open('results.txt','w')

# f.write('Time')
# for i in range(NUMOUTNODES):
	# f.write(' ')
	# f.write('UX_'+str(i))
	# f.write(' ')
	# f.write('UY_'+str(i))
	# f.write(' ')
	# f.write('UZ_'+str(i))

# f.write("\n")

# for i in range (len(time)):
	# f.write(str(time[i]))
	# for k in range(NUMOUTNODES):
		# f.write(' ')
		# f.write(str(UX_CENTERS[k][i]))
		# f.write(' ')
		# f.write(str(UY_CENTERS[k][i]))
		# f.write(' ')
		# f.write(str(UZ_CENTERS[k][i]))
	# f.write("\n")

# f.close()
