from ctypes import byref, cdll, c_int, c_char_p, c_long
from ctypes import POINTER, c_bool, c_double
from numpy import array, float64, zeros

class BATL:

    def __init__(self, libName="/home/abieler/BATL/lib/libWRAPAMR.so"):
        self.lib = cdll.LoadLibrary(libName)

        ################################################################
        # provide argument types of the fortran shared library functions
        self.lib.wrapreadamr_get_nvar.argtypes=[POINTER(c_int)]

        self.lib.wrapreadamr_varnames_.argtypes=[c_char_p, c_long]

        self.lib.wrapreadamr_units_.argtypes=[c_char_p, c_long]

        self.lib.wrapreadamr_read_file_py_.argtypes=[
                                                    POINTER(c_bool),
                                                    POINTER(c_bool),
                                                    c_char_p,
                                                    c_long]

        self.lib.wrapreadamr_domain_limits.argtypes=[POINTER(c_double),
                                                     POINTER(c_double)]

        self.lib.wrapreadamr_get_data_serial.argtypes=[POINTER(c_double),
                                                      POINTER(c_double),
                                                      POINTER(c_bool)]

        self.lib.wrapreadamr_get_more_data.argtypes=[POINTER(c_double),
                                                     POINTER(c_double),
                                                     POINTER(c_bool),
                                                     POINTER(c_int)]
        
        iProc = c_int(0)
        nProc = c_int(1)
        iComm = c_int(91)
        self.lib.wrapreadamr_set_mpi_param(byref(iProc),
                                           byref(nProc),
                                           byref(iComm))

    def load_file(self, dataFileName, isNewGrid=True, isVerbose=False):
        isNewGrid = c_bool(isNewGrid)
        isVerbose = c_bool(isVerbose)
        self.lib.wrapreadamr_read_file_py_(byref(isNewGrid), byref(isVerbose),
                                           dataFileName, len(dataFileName))
        self.nVar = self.get_nVar()
    
    def get_nVar(self):
        nVar = c_int(0)
        self.lib.wrapreadamr_get_nvar(byref(nVar))
        return nVar.value 
    
    def get_data(self, coords):
        coords = array(coords, dtype=float64) 
        isFound = c_bool(False)
        state_V = zeros(self.nVar, dtype=float64) 
        self.lib.wrapreadamr_get_data_serial(coords.ctypes.data_as(POINTER(c_double)),
                                             state_V.ctypes.data_as(POINTER(c_double)),
                                             byref(isFound))
        return state_V, isFound.value

    def get_data_array(self, X):
        '''
        Get interpolated data for multiple locations.
        X is a n x nDim matrix with n = the number of locations
        and nDim (= 1, 2 or 3) the dimension of the model data.
        '''
        n = X.shape[0]
        X = array(X,dtype=float64)
        isFound = c_bool(False)
        state_V = zeros((n, self.nVar), dtype=float64) 

        self.lib.wrapreadamr_get_more_data(X.ctypes.data_as(POINTER(c_double)),
                                           state_V.ctypes.data_as(POINTER(c_double)),
                                           byref(isFound),
                                           byref(c_int(n)))
        return state_V, isFound
    
    def domain_limits(self):
        minCoords = zeros(3,dtype=float64)
        maxCoords = zeros(3,dtype=float64)
        self.lib.wrapreadamr_domain_limits(minCoords.ctypes.data_as(POINTER(c_double)),
                                           maxCoords.ctypes.data_as(POINTER(c_double)))
        return minCoords, maxCoords
    
    def deallocate(self):
        self.lib.wrapreadamr_deallocate()

    def varnames(self):
        varNames = "."*255
        self.lib.wrapreadamr_varnames_(varNames, len(varNames)) 

        return varNames.split()

    def units(self):
        units = "."*255
        self.lib.wrapreadamr_units_(units, len(units))

        return units.split()
