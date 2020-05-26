#!/usr/bin python
from __future__ import print_function
import numpy as np
import read_amr_wrapper

batl = read_amr_wrapper.ReadBATL()
batl.load_file(b"../data/3d__all_3_t00000010_n0000059.idl")
nDim = 3
nPoints = 1000
minDomain, maxDomain = batl.domain_limits()

print("minDomain: ", minDomain)
print("maxDomain: ", maxDomain)

# limit coordinates to be well within the domain boundaries
xMin = minDomain[0] * 0.95
xMax = maxDomain[0] * 0.95

# get number of variables, variable names and units from datafile.
nVars = batl.get_nVar()
print("nVars: ", nVars)

# These work under Python 2, but not python 3
#varNames = batl.varnames()
#for varName in varNames:
#    print(varName)
#units = batl.units()

xyz = np.array([10.0, 10.0, 10.0])
state, isFound = batl.get_data(xyz)
print(state)

X = np.zeros((nPoints, nDim))
xScan = np.linspace(xMin,xMax,nPoints)
X[:,0] = xScan

S, lFound = batl.get_data_array(X)

batl.clean()
