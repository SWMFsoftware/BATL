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
#print(isFound)
print(state)

# Get multiple points from the 3D domain
nDim = 3
nPoints = 10
print("nPoints: ", nPoints)

# Locations (end point is outside)
xMin = minDomain[0] * 0.95
xMax = maxDomain[0] * 1.05
xyz = np.zeros((nPoints, nDim))
xyz[:,0] = np.linspace(xMin,xMax,nPoints)
print("xyz:")
print(xyz)

# Get array of states. Points that are not found are set to zero state
state = batl.get_data_array(xyz)
print("state:")
print(state)

batl.clean()
