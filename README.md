# TreeQSM_Architecture
Functions to extract architectural measures from 3D models of trees based on TLS data

![alt text](https://github.com/TobyDJackson/TreeQSM_Architecture/blob/master/Architecture%20larger%20text.png)

This repo conatins matlab functions to calculate tree architecture models from QSMs. 
If you have TLS point clouds please see https://github.com/InverseTampere/TreeQSM to build QSMs. 

Example QSMs can be found here https://github.com/TobyDJackson/WindAndTrees_FEM/tree/master/QSMs

The script CalculateArchitectures_STRUCT calculates and saves all the architectural measures for QSMs in struct format. The individual functions it calls are contained in Architecture_functions. Other small QSM tools I have found useful are stored in QSM_manipulations. The more interesting architectural measures in this package are: path fraction, crown asymmetry, centre of volume, 
branching order (multiple definitions), mean branching angle, sail area, mass taper exponent. 
