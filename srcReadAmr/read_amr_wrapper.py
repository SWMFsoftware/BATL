from ctypes import byref, cdll, c_int, c_char_p, c_long
from ctypes import POINTER, c_bool, c_double, create_string_buffer
from numpy import array, float64, zeros

class ReadBATL:

    def __init__(self, libName="../lib/libWRAPAMR.so"):
        self.lib = cdll.LoadLibrary(libName)

        ################################################################
        # provide argument types of the fortran shared library functions
        self.lib.wrapamr_get_nvar.argtypes=[POINTER(c_int)]
        self.lib.wrapamr_get_ndim.argtypes=[POINTER(c_int)]
        self.lib.wrapamr_get_namevar.argtypes=[c_char_p, POINTER(c_int)]
        self.lib.wrapamr_get_nameunit.argtypes=[c_char_p, c_int]

        self.lib.wrapamr_read_file.argtypes=[c_char_p,
                                             c_int,
                                             c_int,
                                             c_int]

        self.lib.wrapamr_get_domain.argtypes=[POINTER(c_double),
                                              POINTER(c_double)]

        self.lib.wrapamr_get_data_serial.argtypes=[POINTER(c_double),
                                                      POINTER(c_double),
                                                      POINTER(c_bool)]

        self.lib.wrapamr_get_array_serial.argtypes=[c_int,
                                                    POINTER(c_double),
                                                    POINTER(c_double),
                                                    POINTER(c_bool)]

        self.lib.wrapamr_init_mpi()

    def load_file(self, dataFileName, isNewGrid=True, isVerbose=True):
        isNewGrid = c_int(isNewGrid)
        isVerbose = c_int(isVerbose)
        l = c_int(len(dataFileName))
        self.lib.wrapamr_read_file(dataFileName, l, isNewGrid, isVerbose)
        self.nVar = self.get_nVar()

    def domain_limits(self):
        minCoords = zeros(3,dtype=float64)
        maxCoords = zeros(3, dtype=float64)
        self.lib.wrapamr_get_domain(minCoords.ctypes.data_as(POINTER(c_double)),
                                    maxCoords.ctypes.data_as(POINTER(c_double)))
        return minCoords, maxCoords

    def get_nVar(self):
        nVar = c_int(0)
        self.lib.wrapamr_get_nvar(byref(nVar))
        return nVar.value

    def get_data(self, coords):
        coords = array(coords, dtype=float64)
        isFound = c_bool(False)
        state_V = zeros(self.nVar, dtype=float64)
        self.lib.wrapamr_get_data_serial(coords.ctypes.data_as(POINTER(c_double)),
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

        self.lib.wrapamr_get_array_serial(c_int(n),
                                          X.ctypes.data_as(POINTER(c_double)),
                                          state_V.ctypes.data_as(POINTER(c_double)),
                                          byref(isFound))
        return state_V, isFound

    def clean(self):
        self.lib.wrapamr_clean()

    def varnames(self):
        varNames = "."*500
        l = c_int(0)
        self.lib.wrapamr_get_namevar(varNames, byref(l))
        varNames = varNames[0:l.value]
        return varNames.split()

    def units(self):
        units = "."*255
        self.lib.wrapamr_units_(units, len(units))

        return units.split()
