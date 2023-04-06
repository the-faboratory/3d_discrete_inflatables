
'''
============================================================================
Python code for creating a cylindrical membrane with patches
NB: The units are N and mm (MPa)
Authors: Benjamin Gorissen
Last update: Dec, 2020
Harvard University
email: bgorissen@seas.harvard.edu

Tested on Abaqus 2019
============================================================================
'''

'''
------------------------------------------------------------------------
Import Python-related libraries
------------------------------------------------------------------------
'''

import numpy as np;
from math import *
import copy


'''
------------------------------------------------------------------------
Import Abaqus-related libraries
------------------------------------------------------------------------
'''
from abaqus import *;
from odbAccess import *;
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
import mesh
from time import sleep
from sys import path


# Clear the model database
Mdb()

# Set the directory to the current file directory
os.chdir(os.path.dirname(os.path.realpath('__file__')))
sys.path.append(os.path.dirname(os.path.realpath('__file__')))


session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)


'''
------------------------------------------------------------------------
Helper Functions
------------------------------------------------------------------------
'''

execfile('HelperFunctions.py')

'''
------------------------------------------------------------------------
Simulation options
------------------------------------------------------------------------
'''


#Geometrical parameters
CYL_RADIUS      = 3.18 # 3. Radius of cylinder in mm
CYL_LENGTH      = 292.  # Length of cylinder in mm
CYL_THICKNESS   = 0.183  # 0.33 Cylinder thickness in mm
PATCH_THICKNESS = 0.49 # 0.49
CYL_CIRCUM      = 2*pi*CYL_RADIUS
MESHSIZE        = CYL_CIRCUM/10. #Size of the elements
TOLERANCE       = MESHSIZE/4.

# NUMOUTNODES   = int(ceil(CYL_LENGTH/(2.*MESHSIZE))+1) #Number of central nodes to determine output shape ()
NUMOUTNODES   = 64 #Number of central nodes to determine output shape ()

#------------------------------------------------------------------------
#Patch parameters
#A patch is defined as [[[x1, y1],[x2 y2],[x3 y3], etc.. ], fiber_angle]
#the first list contains the vertices in the xy-plane, 
#where x-coordinates go from 0 to CYL_CIRCUM and y-coordinates go from 0 to CYL_CYL_LENGTH
#width and height are the planar patch dimensions [mm]
#fiber_angle is the fiber angle relative to ... [rad]
#------------------------------------------------------------------------
#To define a rectangular patch you can use the squarepatch command:
#squarepatch(width, height, r,theta,y, patch_angle)
#width and height are the patch dimensions 
#r,theta,y define its position, where r needs to be CYL_RADIUS
#patch_angle is defined as the angle between the longitudinal axis and the height axis of the patch

PATCHES = []

#------------------------------------------------------------------------

PATCHES.append([squarepatch(1.5, 137.2814, CYL_RADIUS, -8.4824/CYL_CIRCUM*(2*pi), 151.9878, 3.6258/180.*pi), 0])



#Selection options
MATERIAL_OPTION = 1   # 1 = Arruda Boyce, 2 = Gent
STEP_OPTION = 3		  # 0=static, 1 = Riks, 2=dynimpl, 3 = dynexpl
LOADING_OPTION = 0	  # 0=pressure, 1=volume
PROLONGED_HOLDING = 1 # 1 = do a second step with reduced mass scaling to let it settle. 

#Material parameters
# AB_LAM           = 2.5    #Arruda boyce labda factor
# AB_MU            = 0.317  #Arruda boyce Shear Modulus
AB_LAM           = 2.3   #Arruda boyce labda factor
AB_MU            = 0.317  #Arruda boyce Shear Modulus
AB_D             = 0.     #Arruda boyce D-factor

AB_LAM_patch     = 2.3

GENT_MU          = 0.3
GENT_K           = 2000. 
GENT_JM			 = 55.


# Neo - Hookean Material
NH_C1_Patch = 0.65 
NH_D1 = 0.

STIFFNESS_FACTOR = 1.1    #How much stiffer the patch is compared to the balloon (multiplier on MU) 

RHO        = 1e-9		  #Density
RHO_Patch  = 1e-8		  #Density
DAMP_ALPHA = 10.		  #Reighleigh Damping
DAMP_BETA  = 0. 		  #Reighleigh Damping

#Loading parameters
STEPNAME_VOLUME ='Inflation-Step'
DT = 3 # The change in temperature for the inflation of the actuator

STEPNAME_PRESSURE ='Pressure-Step'
PRESSURE = 0.02; #Pressure loading in MPa

# Air properties (no real effect on results)
SPECIFIC_HEAT = 29.2
GASCONSTANT = 8.314
MOLWEIGHT = 2.897e-05; #Molecular weight of gas in Mg/mol

if LOADING_OPTION==0: # 0=pressure, 1=volume
	STEPNAME = STEPNAME_PRESSURE
else:
	STEPNAME = STEPNAME_VOLUME
STEPNAME2 = 'Holding-Step'

#Parameters for Riks method
RIKS_INITIAL_ARC = 0.01;
RIKS_MAX_ARC     = 0.1;
RIKS_TOTAL_ARC   = 1.0;
RIKS_MAX_INC     = 400;
RIKS_MIN_ARC     = 1e-06;

#Parameters for Dynamic Implicit
DYNIMP_TIMEPERIOD  = 10.0; 
DYNIMP_MAX_INC     = 400;
DYNIMP_INITIAL_INC = 0.01;
DYNIMP_MIN_INC     =  1e-06; 

#Parameters for Dynamic Explicit
MASSSCALING_FACTOR1 = 100.;
DYNEXPL_TIMEPERIOD1 = 10.;
MASSSCALING_FACTOR2 = 100.;
DYNEXPL_TIMEPERIOD2 = 10.;


#---------------------------------------------------------------------------------------
#Project info
PROJECT_NAME = 'Cylinder';

MODEL_NAME    = PROJECT_NAME+'_model';
PART_NAME     = PROJECT_NAME+'_part';
ASSEMBLY_NAME = PROJECT_NAME+'_assembly';
SKETCH_NAME   = PROJECT_NAME+'_sketch';
CAE_FILE      = PROJECT_NAME+'.cae';
INSTANCE_NAME = PROJECT_NAME+'_instance';
JOB_NAME      = PROJECT_NAME+'_job'

try:
    os.remove(JOB_NAME+'.lck');
    os.remove(JOB_NAME+'.odb');
except OSError:
    pass


'''
------------------------------------------------------------------------
Pre-processing
------------------------------------------------------------------------
'''

#Do some preprocessing on the patches list
NEWPATCHES = []

for PATCH in PATCHES:
	if is_patch_inbounds(PATCH,CYL_CIRCUM,CYL_LENGTH):
		NEWPATCHES.append(PATCH)
	else:
		#NEWPATCHES.append(divide_patch(PATCH,CYL_RADIUS))
		for PATCH2 in divide_patch(PATCH,CYL_RADIUS):
			PATCH2 = make_ccw(PATCH2)
			NEWPATCHES.append(PATCH2)
PATCHES = NEWPATCHES


#Model
Model = mdb.Model(name=MODEL_NAME);

try:
    del mdb.models['Model-1'];
except :
    pass

Model.setValues(absoluteZero=0, universalGas=GASCONSTANT)

#------------------------------------------------------------------------
# PART CYLINDER

#Create flat plate with size of cylinder shell
sketch = Model.ConstrainedSketch(name='__profile__', sheetSize=max(CYL_LENGTH, CYL_CIRCUM))
sketch.rectangle(point1=(0, 0), point2=(CYL_CIRCUM, CYL_LENGTH))
part = Model.Part(name=PART_NAME, dimensionality=THREE_D, type=DEFORMABLE_BODY)
part.BaseShell(sketch=sketch)
del sketch

#Create Sets ### need to change to make it compatible with crossing patches
part.Set(edges=part.edges.findAt(((CYL_CIRCUM/2, 0.0, 0.0), )), name='Set-Bottom')
part.Set(edges=part.edges.findAt(((CYL_CIRCUM/2, CYL_LENGTH, 0.0), )), name='Set-Top')
part.Set(edges=part.edges.findAt(((0.0, CYL_LENGTH/2, 0.0), )), name='Set-Left')
part.Set(edges=part.edges.findAt(((CYL_CIRCUM, CYL_LENGTH/2, 0.0), )), name='Set-Right')

#Partition to get patches
transform = part.MakeSketchTransform(
	sketchPlane=part.faces.findAt(coordinates=(CYL_CIRCUM/2., CYL_LENGTH/2., 0.0), normal=(0.0, 0.0, 1.0)), 
	sketchUpEdge=part.edges.findAt(coordinates=(CYL_CIRCUM, CYL_LENGTH/2., 0.0)), 
	sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))


for PATCH in PATCHES:
	#PATCH is a list that contains a list of coordinates in [0] and the fiber orientation in [1]
	COORD = PATCH[0]
	COORD_CENTR = map(mean, zip(*COORD))

	sketch = Model.ConstrainedSketch(name='__profile__', sheetSize=max(CYL_LENGTH, CYL_CIRCUM), gridSpacing=MESHSIZE, transform=transform)
	sketch.setPrimaryObject(option=SUPERIMPOSE)
	part.projectReferencesOntoSketch(sketch=sketch, filter=COPLANAR_EDGES)

	for index in range(len(COORD)):
		sketch.Line(point1=(COORD[index-1][0], COORD[index-1][1]), point2=(COORD[index][0], COORD[index][1]))

	part.PartitionFaceBySketch(
		sketchUpEdge=part.edges.findAt(coordinates=(CYL_CIRCUM, CYL_LENGTH/2., 0.0)), 
		faces=part.faces.findAt(((COORD_CENTR[0], COORD_CENTR[1], 0.0), )), 
		sketch=sketch)

	del sketch

set_new = part.Set(faces=part.faces.getByBoundingBox(-1.,-1.,-1.,CYL_CIRCUM+1.,CYL_LENGTH+1.,1.), name='Set-Allfaces')

index = 0;
set_patches = []
set_patches.append(set_new)

for PATCH in PATCHES:
	COORD = PATCH[0]
	COORD_CENTR = map(mean, zip(*COORD))

	set_new = part.Set(faces=part.faces.findAt(((COORD_CENTR[0], COORD_CENTR[1], 0.0), )), name='Set-Patch-'+str(index))
	set_patches.append(set_new)
	index = index+1

part.SetByBoolean(name='Set-naked', sets=tuple(set_patches), operation=DIFFERENCE)

#Create Surfaces for pressure
side1Faces = part.faces.getByBoundingBox(-1.,-1.,-1.,CYL_CIRCUM+1.,CYL_LENGTH+1.,1.)
part.Surface(side1Faces=side1Faces, name='Surf-Cylinder')



#Create material
if MATERIAL_OPTION ==1: #Arruda Boyce
	Model.Material(name='FlexibleMaterial')
	Model.materials['FlexibleMaterial'].Hyperelastic(
        materialType=ISOTROPIC, testData=OFF, type=ARRUDA_BOYCE, 
        volumetricResponse=VOLUMETRIC_DATA, table=((AB_MU, AB_LAM, AB_D), ))


	Model.Material(name='RigidMaterial')

	Model.materials['RigidMaterial'].Hyperelastic(
		materialType=ISOTROPIC, testData=OFF, type=NEO_HOOKE, 
        volumetricResponse=VOLUMETRIC_DATA, table=((NH_C1_Patch, NH_D1), ))
        
	#Model.materials['RigidMaterial'].Hyperelastic(
#		materialType=ISOTROPIC, testData=OFF, type=ARRUDA_BOYCE, 
 #       volumetricResponse=VOLUMETRIC_DATA, table=((STIFFNESS_FACTOR*AB_MU, AB_LAM_patch, AB_D), ))

if MATERIAL_OPTION ==2: #Gent
	if STEP_OPTION == 3: #dynexpl
		Model.Material(name='FlexibleMaterial')
		Model.materials['FlexibleMaterial'].Hyperelastic(
		    materialType=ISOTROPIC, testData=OFF, type=USER, 
		    moduliTimeScale=INSTANTANEOUS, properties=3, table=((GENT_MU, GENT_K, GENT_JM), ))

		Model.Material(name='RigidMaterial')
		Model.materials['RigidMaterial'].Hyperelastic(
			materialType=ISOTROPIC, testData=OFF, type=USER, 
		    moduliTimeScale=INSTANTANEOUS, properties=3, table=((STIFFNESS_FACTOR*GENT_MU, GENT_K, GENT_JM), ))
	else:
		Model.Material(name='FlexibleMaterial')
		Model.materials['FlexibleMaterial'].Hyperelastic(
		    materialType=ISOTROPIC, testData=OFF, type=USER, 
		    moduliTimeScale=INSTANTANEOUS, properties=3, table=((GENT_MU, GENT_JM), ))

		Model.Material(name='RigidMaterial')
		Model.materials['RigidMaterial'].Hyperelastic(
			materialType=ISOTROPIC, testData=OFF, type=USER, 
		    moduliTimeScale=INSTANTANEOUS, properties=3, table=((STIFFNESS_FACTOR*GENT_MU, GENT_JM), ))


Model.materials['RigidMaterial'].Density(table=((RHO, ), ))
Model.materials['FlexibleMaterial'].Density(table=((RHO_Patch, ), ))
Model.materials['RigidMaterial'].Damping(alpha=DAMP_ALPHA, beta=DAMP_BETA)
Model.materials['FlexibleMaterial'].Damping(alpha=DAMP_ALPHA, beta=DAMP_BETA)



#Create and assign sections
Model.MembraneSection(name='Section-Flexible', 
    material='FlexibleMaterial', thicknessType=UNIFORM, thickness=CYL_THICKNESS, 
    thicknessField='', poissonDefinition=DEFAULT)

Model.MembraneSection(name='Section-Rigid', 
    material='RigidMaterial', thicknessType=UNIFORM, thickness=PATCH_THICKNESS, 
    thicknessField='', poissonDefinition=DEFAULT)

part.SectionAssignment(region=part.sets['Set-naked'], sectionName='Section-Flexible', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

index = 0;
for PATCH in PATCHES:
	part.SectionAssignment(region=part.sets['Set-Patch-'+str(index)], sectionName='Section-Rigid', offset=0.0, 
		offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
	index = index+1

#Partition to get outputnodes
for index in range(NUMOUTNODES-2):
	datum_temp = part.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=(index+1)*CYL_LENGTH/(NUMOUTNODES-1))
	pickedFaces = part.faces
	part.PartitionFaceByDatumPlane(datumPlane=part.datums[datum_temp.id], faces=pickedFaces)





#Mesh flat plane
pickedRegions = part.faces.getByBoundingBox(-1.,-1.,-1.,CYL_CIRCUM+1.,CYL_LENGTH+1.,1.)
# part.setMeshControls(elemShape=TRI, regions=pickedRegions)
part.setMeshControls(regions=pickedRegions, elemShape=QUAD_DOMINATED)
part.seedPart(size=MESHSIZE, deviationFactor=0.1, minSizeFactor=0.1)

if STEP_OPTION==3: #Explicit
	elemType1 = mesh.ElemType(elemCode=M3D3, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF)
	elemType2 = mesh.ElemType(elemCode=M3D4, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF)
	part.setElementType(regions=tuple(part.faces), elemTypes=(elemType1, elemType2))
else:
	part.setElementType(
		elemTypes=(
			ElemType(elemCode=S4R, elemLibrary=STANDARD), 
			ElemType(elemCode=S3, elemLibrary=STANDARD)), 
		regions=tuple(part.faces))
part.generateMesh()


#Create set of output planes
for index in range(NUMOUTNODES):
	offset=(index)*CYL_LENGTH/(NUMOUTNODES-1)
	pickedNodes = part.nodes.getByBoundingBox(TOLERANCE,offset-TOLERANCE,-TOLERANCE,CYL_CIRCUM+TOLERANCE, offset+TOLERANCE,TOLERANCE)
	part.Set(nodes=pickedNodes, name='Set-output-'+str(index))



#Deform mesh to cylinder
deformation = lambda x: plane_to_cyl(x,CYL_RADIUS) 
partNodes = part.nodes
for partNode in partNodes:
    nodeCoords = partNode.coordinates
    nodeDeformation = deformation(nodeCoords[0])
    # update the coordinate of the node in the mesh
    part.editNode(nodes=partNode, offset1=nodeDeformation[0], offset2=nodeDeformation[1], offset3=nodeDeformation[2])

#------------------------------------------------------------------------
# PART DISK

#Create top disk
sketch2 = Model.ConstrainedSketch(name='__profile2__', sheetSize=CYL_CIRCUM)
coordinates = []
index = 0
for i in range(0,len(part.sets['Set-Top'].nodes)):
	TopCoordinate_1 = part.sets['Set-Top'].nodes[i-1].coordinates
	TopCoordinate_2 = part.sets['Set-Top'].nodes[i].coordinates
	sketch2.Line(point1=(TopCoordinate_1[0], TopCoordinate_1[2]), point2=(TopCoordinate_2[0], TopCoordinate_2[2]))
	coordinates.append(TopCoordinate_1)
	index = index+1
part_disk = Model.Part(name='Disk_part', dimensionality=THREE_D, type=DEFORMABLE_BODY)
part_disk.BaseShell(sketch=sketch2)
del Model.sketches['__profile2__']

# part.setMeshControls(elemShape=TRI, regions=pickedRegions)
part.setMeshControls(regions=pickedRegions, elemShape=QUAD_DOMINATED)

#Mesh top disk
pickedRegions = part_disk.faces
# part_disk.setMeshControls(elemShape=TRI, regions=pickedRegions)
part_disk.setMeshControls(regions=pickedRegions, elemShape=QUAD_DOMINATED)
part_disk.seedPart(size=MESHSIZE, deviationFactor=0.1, minSizeFactor=0.1)
if STEP_OPTION==3: #Explicit
	elemType1 = mesh.ElemType(elemCode=M3D3, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF)
	elemType2 = mesh.ElemType(elemCode=M3D4, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF)
	part_disk.setElementType(regions=tuple(part_disk.faces), elemTypes=(elemType1, elemType2))
else:
	part_disk.setElementType(
		elemTypes=(
			ElemType(elemCode=S4R, elemLibrary=STANDARD), 
			ElemType(elemCode=S3, elemLibrary=STANDARD)), 
		regions=tuple(part_disk.faces))
part_disk.generateMesh()



#Create set and surface on disk 
node_list = part_disk.nodes[0:index]
part_disk.Set(nodes=node_list, name='Set-Top-Disk')
part_disk.Surface(side1Faces=part_disk.faces, name='Surf-Disk')
part_disk.Set(faces=part_disk.faces.findAt(((0.0, 0.0, 0.0), )), name='Set-Disk-All')

#Assign section to disk 
part_disk.SectionAssignment(region=part_disk.sets['Set-Disk-All'], sectionName='Section-Flexible', offset=0.0, 
		offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

#------------------------------------------------------------------------
# ASSEMBLY

#Assembly
Model.rootAssembly.DatumCsysByDefault(CARTESIAN)
assembly = Model.rootAssembly
instance = assembly.Instance(dependent=ON, name=ASSEMBLY_NAME, part=part)
instance_disk = assembly.Instance(dependent=ON, name='Disk_assembly', part=part_disk)
assembly.regenerate()

assembly.rotate(instanceList=('Disk_assembly', ), axisPoint=(0.0, 0.0, 0.0), 
        axisDirection=(1.0, 0.0, 0.0), angle=90.0)
assembly.translate(instanceList=('Disk_assembly', ), vector=(0.0, CYL_LENGTH, 0.0))

#Create Set of everything on cylinder
assembly.Set(faces=instance.faces, name='Set-AllCylinder')

#Tie Nodes together to seal cylinder to itself
ranNodesRight=range(0,len(instance.sets['Set-Right'].nodes))
ind1=1
for i in range(0,len(instance.sets['Set-Left'].nodes)):
	LeftCoordinate = instance.sets['Set-Left'].nodes[i].coordinates
	assembly.Set(name='Node-Left-'+str(ind1), nodes=instance.sets['Set-Left'].nodes[i:i+1])
	for j in ranNodesRight:
		RightCoordinate = instance.sets['Set-Right'].nodes[j].coordinates
		if ((sqrt((LeftCoordinate[0]-RightCoordinate[0])**2)) < TOLERANCE and sqrt((LeftCoordinate[1]-RightCoordinate[1])**2) < TOLERANCE and sqrt((LeftCoordinate[2]-RightCoordinate[2])**2) < TOLERANCE): #fabs(sqrt((LeftCoordinate[0]-RightCoordinate[0])**2+(LeftCoordinate[1]-RightCoordinate[1])**2+(LeftCoordinate[2]-RightCoordinate[2])**2)-LatticeVec[0][0]) < TOL and 
			assembly.Set(name='Node-Right-'+str(ind1), nodes=instance.sets['Set-Right'].nodes[j:j+1])
			Model.Equation(name='TIE-X'+'-'+str(ind1),
				terms=((1.0,'Node-Left-'+str(ind1), 1),(-1.0, 'Node-Right-'+str(ind1), 1)))
			Model.Equation(name='TIE-Y'+'-'+str(ind1),
				terms=((1.0,'Node-Left-'+str(ind1), 2),(-1.0, 'Node-Right-'+str(ind1), 2)))
			Model.Equation(name='TIE-Z'+'-'+str(ind1),
				terms=((1.0,'Node-Left-'+str(ind1), 3),(-1.0, 'Node-Right-'+str(ind1), 3)))

			ind1=ind1+1
			break

#Tie Nodes together to seal cylinder and disk
ranNodesTopCyl=range(0,len(instance.sets['Set-Top'].nodes))
ind1=1
for i in range(0,len(instance_disk.sets['Set-Top-Disk'].nodes)):
	DiskCoordinate = instance_disk.sets['Set-Top-Disk'].nodes[i].coordinates
	assembly.Set(name='Node-TopDisk-'+str(ind1), nodes=instance_disk.sets['Set-Top-Disk'].nodes[i:i+1])
	for j in ranNodesTopCyl:
		CylCoordinate = instance.sets['Set-Top'].nodes[j].coordinates
		if ((sqrt((DiskCoordinate[0]-CylCoordinate[0])**2)) < TOLERANCE and sqrt((DiskCoordinate[1]-CylCoordinate[1])**2) < TOLERANCE and sqrt((DiskCoordinate[2]-CylCoordinate[2])**2) < TOLERANCE): 
			assembly.Set(name='Node-TopCyl-'+str(ind1), nodes=instance.sets['Set-Top'].nodes[j:j+1])
			Model.Equation(name='TOPTIE-X'+'-'+str(ind1),
				terms=((1.0,'Node-TopDisk-'+str(ind1), 1),(-1.0, 'Node-TopCyl-'+str(ind1), 1)))
			Model.Equation(name='TOPTIE-Y'+'-'+str(ind1),
				terms=((1.0,'Node-TopDisk-'+str(ind1), 2),(-1.0, 'Node-TopCyl-'+str(ind1), 2)))
			Model.Equation(name='TOPTIE-Z'+'-'+str(ind1),
				terms=((1.0,'Node-TopDisk-'+str(ind1), 3),(-1.0, 'Node-TopCyl-'+str(ind1), 3)))

			ind1=ind1+1
			break


#Create fluid Cavity
assembly.ReferencePoint(point=(0.0, 0.0, 0.0))
assembly.features.changeKey(fromName='RP-1', toName='CavityRP')
cavityPointId = assembly.features['CavityRP'].id
assembly.Set(name='CavityRP', referencePoints=(Model.rootAssembly.referencePoints[cavityPointId], ))


if LOADING_OPTION==0: #0=pressure, 1=volume
	if STEP_OPTION==3:
		Model.FluidCavityProperty(name='PneumaticCavityProp', 
	        definition=PNEUMATIC, molecularWeight=MOLWEIGHT,
	        useCapacity=True, capacityTable=((SPECIFIC_HEAT, 0.0, 0.0, 0.0, 0.0), ))
	else:
		Model.FluidCavityProperty(name='PneumaticCavityProp', 
	        definition=PNEUMATIC, molecularWeight=MOLWEIGHT,
	        useCapacity=False)
elif LOADING_OPTION==1:
	if STEP_OPTION==3: #Explicit
		Model.FluidCavityProperty(
			expansionTable=((1.0, ), ), 
			fluidDensity=1e-09, 
			name='PneumaticCavityProp', 
			useExpansion=True,
			useBulkModulus=True, 
			bulkModulusTable=((2150.0, ),))
	else:
		Model.FluidCavityProperty(
			expansionTable=((1.0, ), ), 
			fluidDensity=1e-09, 
			name='PneumaticCavityProp', 
			useExpansion=True)
else:
	print('Error in defining step')

assembly.SurfaceByBoolean(name='Surf-Cavity', surfaces=(instance.surfaces['Surf-Cylinder'], instance_disk.surfaces['Surf-Disk'], ))


Model.FluidCavity(name='PneumaticCavity', 
    createStepName='Initial', cavityPoint=assembly.sets['CavityRP'], cavitySurface=assembly.surfaces['Surf-Cavity'], 
    interactionProperty='PneumaticCavityProp')

Model.VelocityBC(name='BC-Velocity', createStepName='Initial', 
        region=instance.sets['Set-Bottom'], v1=0.0, v2=0.0, v3=0.0, vr1=UNSET, vr2=UNSET, vr3=UNSET, 
        amplitude=UNSET, localCsys=None, distributionType=UNIFORM, 
        fieldName='')


#Create Step

if STEP_OPTION == 0: #0=static, 
	Model.StaticStep(name=STEPNAME, previous='Initial', 
    maxNumInc=100, initialInc=0.0001, minInc=1e-07, maxInc=0.1, nlgeom=ON)
elif STEP_OPTION == 1: #1 = Riks, 
	Model.StaticRiksStep(name=STEPNAME, previous='Initial', 
		maxNumInc=RIKS_MAX_INC, initialArcInc=RIKS_INITIAL_ARC, minArcInc=RIKS_MIN_ARC, maxArcInc=RIKS_MAX_INC, 
		totalArcLength=RIKS_TOTAL_ARC, nlgeom=ON)
elif STEP_OPTION == 2: #2=dynimpl
	Model.ImplicitDynamicsStep(name=STEPNAME, 
        previous='Initial', timePeriod=DYNIMP_TIMEPERIOD, maxNumInc=DYNIMP_MAX_INC, 
        application=DEFAULT, initialInc=DYNIMP_INITIAL_INC, minInc=DYNIMP_MIN_INC, nohaf=OFF, 
        amplitude=RAMP, alpha=DEFAULT, initialConditions=DEFAULT, nlgeom=ON)
elif STEP_OPTION == 3: #3=dynamic explicit
	Model.ExplicitDynamicsStep(name=STEPNAME, 
        previous='Initial', massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, 
        MASSSCALING_FACTOR1, 0.0, None, 0, 0, 0.0, 0.0, 0, None), ), timePeriod=DYNEXPL_TIMEPERIOD1, scaleFactor=1.0, 
        linearBulkViscosity=0.06, quadBulkViscosity=1.2,improvedDtMethod=ON)
	Model.steps[STEPNAME].setValues(massScaling=((
        SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 0.0001, BELOW_MIN, 1000, 
        0, 0.0, 0.0, 0, None), ), improvedDtMethod=ON)
	if PROLONGED_HOLDING == 1:
		Model.ExplicitDynamicsStep(name=STEPNAME2, 
        	previous=STEPNAME, massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, 
        	MASSSCALING_FACTOR2, 0.0, None, 0, 0, 0.0, 0.0, 0, None), ), timePeriod=DYNEXPL_TIMEPERIOD2, scaleFactor=1.0, 
        	linearBulkViscosity=0.06, quadBulkViscosity=1.2,improvedDtMethod=ON)
        Model.steps[STEPNAME2].setValues(
        	timeIncrementationMethod=FIXED_USER_DEFINED_INC, userDefinedInc=0.0001, 
        	massScaling=((DISABLE_THROUGHOUT_STEP, MODEL, THROUGHOUT_STEP, 0.0, 
        	0.0, None, 1, 0, 0.0, 0.0, 0, None), (SEMI_AUTOMATIC, MODEL, 
        	AT_BEGINNING, MASSSCALING_FACTOR2, 0.0, None, 0, 0, 0.0, 0.0, 0, None)), 
        	improvedDtMethod=ON)

else:
	print('Error in defining step')

# Creating load
if STEP_OPTION==0 or STEP_OPTION==1: # 0=static, 1 = Riks, 2=dynimpl
	if LOADING_OPTION==0: #0=pressure, 1=volume
		Model.Pressure(createStepName=STEPNAME, distributionType=UNIFORM, magnitude=PRESSURE, name='pload', region=assembly.surfaces['Surf-Cavity'])
	else:
		Model.Temperature(name='Initial_Temp', 
			createStepName='Initial', region=assembly.sets['CavityRP'], distributionType=UNIFORM, 
			crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(0.0, ))
		Model.Temperature(name='Final_Temp', 
			createStepName=STEPNAME_VOLUME, region=assembly.sets['CavityRP'], 
			distributionType=UNIFORM, 
			crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(DT, ))
if STEP_OPTION==2:
	if LOADING_OPTION==0:
		Model.Pressure(createStepName=STEPNAME, distributionType=UNIFORM, magnitude=PRESSURE, name='pload', region=assembly.surfaces['Surf-Cavity'])
	else:
		Model.Temperature(name='Initial_Temp', 
			createStepName='Initial', region=assembly.sets['CavityRP'], distributionType=UNIFORM, 
			crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(0.0, ))
		Model.SmoothStepAmplitude(name='Amplitude_inflation', timeSpan=STEP, data=((0.0, 0.0), (DYNIMP_TIMEPERIOD, 1.0)))
		Model.Temperature(name='Final_Temp', 
			createStepName=STEPNAME_VOLUME, region=assembly.sets['CavityRP'], 
			distributionType=UNIFORM, 
			crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, 
			magnitudes=(DT, ), amplitude='Amplitude_inflation')
if STEP_OPTION==3:
	if LOADING_OPTION==0:
		# for index in range(NUMOUTNODES-1):
		# 	Model.SmoothStepAmplitude(name='Amplitude_inflation'+str(index), timeSpan=STEP, data=((0.0, 0.0), (2.*float(index)+0.01, 0.0), (DYNEXPL_TIMEPERIOD1+2.*float(index), 1.0), (2*DYNEXPL_TIMEPERIOD1, 1.0)))
		# 	Model.Pressure(createStepName=STEPNAME, distributionType=UNIFORM, magnitude=PRESSURE, amplitude='Amplitude_inflation'+str(index), name='pload'+str(index), region=instance.surfaces['Surf-Cylinder'+str(index)])
		Model.SmoothStepAmplitude(name='Amplitude_inflation', timeSpan=STEP, data=((0.0, 0.0), (DYNEXPL_TIMEPERIOD1-2., 1.0), (DYNEXPL_TIMEPERIOD1, 1.0)))
		Model.Pressure(createStepName=STEPNAME, distributionType=UNIFORM, magnitude=PRESSURE, amplitude='Amplitude_inflation', name='pload', region=assembly.surfaces['Surf-Cavity'])
		if PROLONGED_HOLDING == 1:
			Model.loads['pload'].deactivate(STEPNAME2)
			Model.Pressure(createStepName=STEPNAME2, distributionType=UNIFORM, magnitude=PRESSURE, amplitude=UNSET, name='pload2', region=assembly.surfaces['Surf-Cavity'])

	else:
		Model.Temperature(name='Initial_Temp', 
			createStepName='Initial', region=assembly.sets['CavityRP'], distributionType=UNIFORM, 
			crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(0.0, ))
		Model.SmoothStepAmplitude(name='Amplitude_inflation', timeSpan=STEP, data=((0.0, 0.0), (DYNEXPL_TIMEPERIOD1, 1.0)))
		Model.Temperature(name='Final_Temp', 
			createStepName=STEPNAME_VOLUME, region=assembly.sets['CavityRP'], 
			distributionType=UNIFORM, 
			crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, 
			magnitudes=(DT, ), amplitude='Amplitude_inflation')
else:
	print('Error in creating loading')


#------------------------------------------------------------------------
# CREATING OUTPUTS

# History Output requests
try:
    del Model.historyOutputRequests['H-Output-1']
except:
    pass

# Model.HistoryOutputRequest(createStepName=STEPNAME, name='PV', region=assembly.sets['CavityRP'], variables=('PCAV', 'CVOL'), numIntervals=2)
Model.HistoryOutputRequest(name='ENERGY_OUTPUT', createStepName=STEPNAME, variables=('ALLAE', 'ALLCD', 'ALLDC', 'ALLDMD', 'ALLFD', 'ALLIE', 'ALLKE', 'ALLPD', 'ALLSE', 'ALLVD', 'ALLWK', 'ALLCW', 'ALLMW', 'ALLPW', 'ETOTAL'))
# Model.HistoryOutputRequest(name='DISPLACEMENTS', 
# 	createStepName=STEPNAME, variables=('U1', 'U2', 'U3'), 
# 	numIntervals=1, region=assembly.sets['Set-AllCylinder'])

if PROLONGED_HOLDING == 0:
	for index in range(NUMOUTNODES):
		Model.HistoryOutputRequest(name='DISPLACEMENT_OUTPUT_'+str(index), 
	        createStepName=STEPNAME, variables=('U1', 'U2', 'U3'), numIntervals=2,
	        region=assembly.sets[ASSEMBLY_NAME+'.Set-output-'+str(index)])
if PROLONGED_HOLDING == 1:
	for index in range(NUMOUTNODES):
		Model.HistoryOutputRequest(name='DISPLACEMENT_OUTPUT_'+str(index), 
	        createStepName=STEPNAME2, variables=('U1', 'U2', 'U3'), numIntervals=2,
	        region=assembly.sets[ASSEMBLY_NAME+'.Set-output-'+str(index)])


#Job
if MATERIAL_OPTION == 1: # 1 = Arruda Boyce, 2 = Gent
	j = mdb.Job(name=JOB_NAME,model=MODEL_NAME,numCpus=4,numDomains=4,multiprocessingMode=DEFAULT,description='cylinder_simulation');
if MATERIAL_OPTION == 2: 
	if STEP_OPTION==3:
		j = mdb.Job(name=JOB_NAME,model=MODEL_NAME,numCpus=4,numDomains=4,multiprocessingMode=DEFAULT,description='cylinder_simulation', 
			userSubroutine='C:\\Users\\u0063816\\OneDrive - Harvard University\\GorissenBen_Collaborations\\15.RebeccaKramer\\3.MaterialFitting\\vumat.for');
	else:
		j = mdb.Job(name=JOB_NAME,model=MODEL_NAME,numCpus=4,numDomains=4,multiprocessingMode=DEFAULT,description='cylinder_simulation', 
			userSubroutine='C:\\Users\\u0063816\\OneDrive - Harvard University\\GorissenBen_Collaborations\\15.RebeccaKramer\\3.MaterialFitting\\uhyper_gent.for');
Model.rootAssembly.regenerate();
mdb.jobs[JOB_NAME].submit(consistencyChecking=OFF)
mdb.jobs[JOB_NAME].waitForCompletion()