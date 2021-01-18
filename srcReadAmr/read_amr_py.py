#!/usr/bin python
from __future__ import print_function
import numpy as np
import read_amr_wrapper

# Initialize the wrapper
# If necessary, pass the library location here. Default is ../lib/libWRAPAMR.so
#batl = read_amr_wrapper.ReadBATL(libName="full_path_to/libWRAPAMR.so")
batl = read_amr_wrapper.ReadBATL()

# Load data file
batl.load_file(b"../data/3d__all_3_t00000010_n0000059.idl")

# Get domain limits
minDomain, maxDomain = batl.domain_limits()

print("minDomain: ", minDomain)
print("maxDomain: ", maxDomain)

# limit coordinates to be well within the domain boundaries
xMin = minDomain[0] * 0.95
xMax = maxDomain[0] * 0.95

# get number of variables, variable names and units from datafile.
nVars = batl.get_nVar()
print("nVars: ", nVars)

# These work under Python 2, but not Python 3
#varNames = batl.varnames()
#for varName in varNames:
#    print(varName)
#units = batl.units()

# get state variables at some location
xyz = np.array([10.0, 10.0, 10.0])
state, isFound = batl.get_data(xyz)
print(state)

# Get multiple points
nDim = 3
nPoints = 10

# Locations
X = np.zeros((nPoints, nDim))
xScan = np.linspace(xMin,xMax,nPoints)
X[:,0] = xScan

state, isFound = batl.get_data_array(X)
#print(isFound)
#print(state)

batl.clean()
