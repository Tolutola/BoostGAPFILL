# BoostGAPFILL

beta version of BoostGAPFILL for gap filling analysis of metabolic network reconstructions implemented in MATLAB


Requirements: MATLAB with COBRA toolbox installed, IBM CPLEX or GUROBI solver (both have free fully functional academic licenses) and any version of Python.

Download the latest version of BoostGapFill from https://github.com/Tolutola/BoostGAPFILL

To see a demo, change the MATLAB working directory to BoostGAPFILL and type the following in the command line:
a)	‘BoostGAPFILL_example1’ to see a simple demo run on iAF1260 model
b)	‘BoostGAPFILL_example2’ to see an extended comparison of BoostGapFill and FASTGAPFILL

The optional settings can be changed in the example scripts. All scripts and function files are in the ‘code’ sub folder. The COBRA models are in the ‘data’ sub folder. The universal reactions, metabolites and stoichiometric matrix are stored as variables in the COBRA model structure. 
