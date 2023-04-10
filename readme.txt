This code accompanies the paper "Programming 3D Curves with Discretely Constrained Cylindrical Inflatables"
R Baines, SK Patiballa, B Gorissen, K Bertoldi, R Kramer‐Bottiglio
Advanced Materials, 2300535

-----------------------------------------------------------------------------------------


fe >>> This directory contains scripts for running the finite element analysis nested in the optimization loop. 

Main script: 
-Main_Design_Opt_IDEP.py
(IMPORTANT: run on a python IDE such as sypder or pycharm for the example, otherwise there might be dependency issues. The script will delete other files in the directory! In our default implementation, this is the C:\\Temp folder. Change in the script to whatever you desire. Upon completion, design outputs are printed to the console. Note that these are scaled outputs, and must be un-scaled by the factors at the end of the Main_Design_Opt_IDEP.py script to get quantities with units.)

Abaqus scripting file: 
-Run.py

Helper functions: 
-HelpreFunctions.py 

Data processing: 
-patched_cylinder_readout.py

Example centerline (peano curve): 
-req_centerline_test_scaled.csv 


ksa >>> This directory contains scripts for running kinematic segmentation algorithm. 

Main script (run this script in the matlab IDE for the example): 
-main.m

Helper functions: 
-arclength.m
-cumulative_arc_len.m
-plot_tube.m

Example centerlines: 
-peano_curve.csv 
-pancake.csv 
