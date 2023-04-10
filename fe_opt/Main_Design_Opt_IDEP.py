# Main file for the shape optimization of the inflatables with discrete elastic patches
'''
-----------------------------------------------------------------------------------------------------
This is the main file for design optmization in abaqus
Author : Sree Patiballa, Yale Universtiy
Current Affiation: University of Alabama, spatiballa@ua.edu
References : Add all the necessary references here
-> Analysis code:  Ben G, from Prof. Bertoldi's group at Harvard
-> CMA-ES documentation : http://cma.gforge.inria.fr/cmaes_sourcecode_page.html#practical

USAGE:
This code runs the forward simulation without GUI and also extracts data from .odb file.
Then the data is postprocessed to get the centerline information at chosen nodes.
This information is given to the optimization algorithm to reduce the RMSE (root mean square error)
b/w the analysis and the required shape (curvature and torsion) in an iterative manner.

-----------------------------------------------------------------------------------------------------
'''

from IPython import get_ipython
get_ipython().magic('reset -sf') # resets/clears the variables

'''
-----------------------------------------------------------------------------------------------------
Import necessary modules, To use these libraries makes sure they are installed
-----------------------------------------------------------------------------------------------------
'''
import numpy as np
from numpy import linalg
from numpy import genfromtxt 
import matplotlib.pyplot as plt
import math
import os
import glob
import csv
import math
import time
from scipy.optimize import differential_evolution
import cma
from HelperFunctions import csv2list
from HelperFunctions import frenet


def tic():
    #Homemade version of matlab tic and toc functions
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print("Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds.")
    else:
        print("Toc: start time not set")


COUNTER = 0

'''
-----------------------------------------------------------------------------------------------------
Clean Up the temp folder after each iteration
Change the names of the files or the folder accordingly
CAUTION : Make sure the Temp or any other folder is backed up. This function clears everything in that folder !!!
-----------------------------------------------------------------------------------------------------
'''
def cleanUp():

    cleanupFiles = glob.glob("C:\\Temp\\*.*") #Cleaning previous iteration files
    try:
      cleanupFiles.remove("C:\\Temp\\patched_cylinder_readoutput_sree.py")
    except ValueError:
      pass
    try:
      cleanupFiles.remove("C:\\Temp\\Main_Design_Opt_IDEP.py")
    except ValueError:
      pass
    try:
      cleanupFiles.remove("C:\\Temp\\Run.py")
    except ValueError:
      pass
    try:
      cleanupFiles.remove("C:\\Temp\\Cylinder_job.sta")
    except ValueError:
      pass
    try:
        cleanupFiles.remove("C:\\Temp\\centerline_coord.csv")
    except ValueError:
        pass
    try:
        cleanupFiles.remove("C:\\Temp\\req_centerline_test_scaled.csv")
    except ValueError:
        pass
    try:
        cleanupFiles.remove("C:\\Temp\\HelperFunctions.py")
    except ValueError:
        pass

    os.chdir("C:\\Temp")
    for cleanupFile in cleanupFiles:
        os.remove(cleanupFile)


'''
-----------------------------------------------------------------------------------------------------
This functions opens up the model/analysis file and
modifies the necessary information and writes the design values of the design variables
-----------------------------------------------------------------------------------------------------
'''
def getSimData(x):
    failure_count = 0 # This is used to check the status of the analysis. 0 implies no errors and 1 or more implies errors

    with open('C:\\Temp\\Run.py','r') as f:
        data = f.readlines()
    f.close()

    # this line changes the patch definition in the Run.py file. Edit all the patches accordingly
    data[107]='PATCHES.append([squarepatch(%f, %f, CYL_RADIUS, 4.0/CYL_CIRCUM*(2*pi), %f, 0.0/180.*pi), 0])\n'%(x[0],x[1],x[2])
    data[108]='PATCHES.append([squarepatch(%f, %f, CYL_RADIUS, 4.0/CYL_CIRCUM*(2*pi), %f, 0.0/180.*pi), 0])\n'%(x[3],x[4],x[5])
    data[109]='PATCHES.append([squarepatch(%f, %f, CYL_RADIUS, 8.99/CYL_CIRCUM*(2*pi), %f, 0.0/180.*pi), 0])\n'%(x[6],x[7],x[8])
    data[110]='PATCHES.append([squarepatch(%f, %f, CYL_RADIUS, 8.99/CYL_CIRCUM*(2*pi), %f, 0.0/180.*pi), 0])\n'%(x[9],x[10],x[11])
    data[111]='PATCHES.append([squarepatch(%f, %f, CYL_RADIUS, 4.0/CYL_CIRCUM*(2*pi), %f, 0.0/180.*pi), 0])\n'%(x[12],x[13],x[14])
    data[112]='PATCHES.append([squarepatch(%f, %f, CYL_RADIUS, 4.0/CYL_CIRCUM*(2*pi), %f, 0.0/180.*pi), 0])\n'%(x[15],x[16],x[17])


    with open('C:\\Temp\\Run.py','w') as file:
        file.writelines(data)
    file.close()

    # calling the cleanup function to clean up the temp folder
    cleanUp()

    # run the analysis file in abaqus with no GUI
    os.system("abaqus cae noGUI=Run.py")

    time.sleep(0.5)

    try:
        with open('C:\\Temp\\Cylinder_job.sta','r') as f_check:
            data_check1 = f_check.readlines()
        f_check.close()


        if(data_check1[-2]== '  THE ANALYSIS HAS COMPLETED SUCCESSFULLY\n'):
            print("Success!")
            # need to run a postprocessing script here and write the data to a csv file
            os.system("abaqus cae noGUI=patched_cylinder_readoutput_sree.py")
        else:
            print("ERROR IN THE ANALYSIS FILE :( ")
            failure_count+=1

    except:
        failure_count+=1


    return failure_count


'''
-----------------------------------------------------------------------------------------------------
Definition of the objective function (with constraints as penalties, if needed)
-----------------------------------------------------------------------------------------------------
'''
def objective_func(x):

     global COUNTER

    # scale up the design variabels to give into ABAQUS
    #x_mul = [0.05, 0.02, 0.25, 0.04, 0.1]
    #x = x_mul * x
     x_mul = [1/2, 1/20, 1/50, 1/2, 1/20, 1/50, 1/2, 1/20, 1/100, 1/2, 1/20, 1/100, 1/2, 1/20, 1/200, 1/2, 1/20, 1/200]
     x_mul = np.asarray(x_mul)
     x_unscaled = x / x_mul
     x_unscaled = np.ndarray.tolist(x_unscaled)
     FC = getSimData(x_unscaled) # Failure count


    # change to the folder where readData.py script saved the data
     os.chdir("C:\\Temp")

     f = 100 # start the obj function value with 100

     if(FC==0):

        # TODO : Convert the list to a numpy array
        Analysis_Shape = np.genfromtxt('centerline_coord.csv', delimiter=',')
        
        [k1,t1] = frenet(Analysis_Shape[:,0],Analysis_Shape[:,1],Analysis_Shape[:,2])

         # Give the required shape array here
        Req_Shape = np.genfromtxt('req_centerline_test_scaled.csv', delimiter=',')
        
        [k2,t2] = frenet(Req_Shape[:,0],Req_Shape[:,1],Req_Shape[:,2])

         # objective funciton - RMSE between the shape array obtained from analysis and required shape array - divide it by function value with initial x0
        # f = np.linalg.norm(Req_Shape - Analysis_Shape) / 1050
        f = (np.linalg.norm(k1 - k2) + np.linalg.norm(t1 - t2)) / 0.829

        print(f)


        # constraints are defined on the bounds of the design domain and isotropy constraints
        if(x[0] < 0.5 or x[0] > 1.5) or (x[1] < 0.5 or x[1] > 1.5) or (x[2] < 0.8 or x[2] > 1.0) \
            or (x[3] < 0.5 or x[3] > 1.5) or (x[4] < 0.5 or x[4] > 1.5) or (x[5] < 1.6 or x[5] > 1.8) \
                or (x[6] < 0.5 or x[6] > 1.5) or (x[7] < 0.5 or x[7] > 1.5) or (x[8] < 1.2 or x[8] > 1.3) \
                    or (x[9] < 0.5 or x[9] > 1.5) or (x[10] < 0.5 or x[10] > 1.5) or (x[11] < 1.6 or x[11] > 1.7)\
                        or (x[12] < 0.5 or x[12] > 1.5) or (x[13] < 0.5 or x[13] > 1.5) or (x[14] < 1.0 or x[14] > 1.05) \
                            or (x[15] < 0.5 or x[15] > 1.5) or (x[16] < 0.5 or x[16] > 1.5) or (x[17] < 1.2 or x[17] > 1.25):
             f = 100
             print('Out of Bounds')


     COUNTER = COUNTER + 1
     print(x,f,COUNTER)

     return f





'''
-----------------------------------------------------------------------------------------------------
Main file for the optmization algorithm with bounds on the design variables
-----------------------------------------------------------------------------------------------------
'''
if __name__ == '__main__':

  # bounds = [(0.0,0.5), (0.0,0.5), (0.0,0.5), (0.0,0.5), (0.0,0.5)]

  x0 = [1.5, 17, 43.93, 1.5, 17, 87.88, 1.5, 17, 128.21, 1.5, 17, 163.87, 1.5, 18.46, 204.84, 1.5, 18.44, 248.78]
  x0 = np.asarray(x0)
  x_mul = [1/2, 1/20, 1/50, 1/2, 1/20, 1/50, 1/2, 1/20, 1/100, 1/2, 1/20, 1/100, 1/2, 1/20, 1/200, 1/2, 1/20, 1/200]
  x_mul = np.asarray(x_mul)
  x0_scaled = x_mul * x0
  x0_scaled = np.ndarray.tolist(x0_scaled) # [0.75  , 1.    , 0.89  , 0.75  , 1.    , 1.702 , 0.75  , 1.    ,   1.257 , 0.75  , 1.    , 1.663 , 0.75  , 1.    , 1.0345, 0.75  , 1.    , 1.2375]

  tic()
  res = cma.fmin(objective_func, x0_scaled, 0.01, options={'verb_disp':1, 'maxiter': 5})
  toc()

  print(res[0])




