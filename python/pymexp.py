
import ctypes
import ctypes.util

mexppath = ctypes.util.find_library("mexp")
if mexppath:
    # Found, load it:
    libmexp = ctypes.CDLL(mexppath)
else:
    # Not found, try loading Windows DLLs directly:
    try:
        # Try loading 32bit library
        libmexp = ctypes.CDLL("mexp32.dll")
    except OSError:
        # Didn't load. How about the 64bit version:
        libmexp = ctypes.CDLL("mexp64.dll")
        # Throws another OSError if fails, I don't care at this point


_work_alloc = libmexp.mexp_work_alloc
_work_alloc.argtypes = []
_work_alloc.restype = ctypes.c_void_p

_work_free = libmexp.mexp_work_free
_work_free.argtypes = [ctypes.c_void_p]
_work_free.restype = None

_map_coefs = libmexp.mexp_map_coefs
_map_coefs.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double,
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.c_void_p, ctypes.c_void_p]
_map_coefs.restype = ctypes.c_int

_two_objectives = libmexp.mexp_two_objectives
_two_objectives.argtypes = [ctypes.c_double, ctypes.c_double, 
                            ctypes.POINTER(ctypes.c_double),
                            ctypes.POINTER(ctypes.c_double),
                            ctypes.POINTER(ctypes.c_double),
                            ctypes.POINTER(ctypes.c_double),
                            ctypes.c_void_p, ctypes.c_void_p]

     
class NumericException(Exception):
    pass
                 
class Pars(ctypes.Structure):
    _fields_ = [
        ("phi1", ctypes.c_double),
        ("phi2", ctypes.c_double),
        ("varphi", ctypes.c_double),
        ("k0_1", ctypes.c_double),
        ("k0_2", ctypes.c_double),
        ("L0_A", ctypes.c_double),
        ("L0_B", ctypes.c_double),
        ("rho", ctypes.c_double),
        ("xi", ctypes.c_double),
        ("mu", ctypes.c_double),
        ("zeta", ctypes.c_double),
        ("psi", ctypes.c_double),
        ("chi", ctypes.c_double),
        ("gamma", ctypes.c_double),
        ("q0", ctypes.c_double),
        ("q1", ctypes.c_double),
        ("X1", ctypes.c_double)
    ]

class Mexp:
    def __init__(self, *args, **kwargs):
        self._work = _work_alloc()
        self.pars = Pars(*args, **kwargs)

    def __del__(self):
        _work_free(self._work)

    def map_coefs(self, eps_imp_1, eps_fcl_1, eps_imp_2):
        c_eps_imp_1 = ctypes.c_double(eps_imp_1)
        c_eps_fcl_1 = ctypes.c_double(eps_fcl_1)
        c_eps_imp_2 = ctypes.c_double(eps_imp_2)
        c_Lambda_1 = ctypes.c_double()
        c_Kappa_1 = ctypes.c_double()
        c_Iota_1 = ctypes.c_double()
        c_Lambda_2 = ctypes.c_double()
        status = _map_coefs(c_eps_imp_1, c_eps_fcl_1, c_eps_imp_2,
                            ctypes.byref(c_Iota_1),
                            ctypes.byref(c_Kappa_1),
                            ctypes.byref(c_Lambda_1),
                            ctypes.byref(c_Lambda_2),
                            ctypes.byref(self.pars),
                            self._work)
        if status != 0:
            raise NumericException()

        return (c_Iota_1.value, c_Kappa_1.value, c_Lambda_1.value, c_Lambda_2.value)

    def objectives(self, k_fcl_A, k_fcl_B):
        c_k_fcl_A = ctypes.c_double(k_fcl_A)
        c_k_fcl_B = ctypes.c_double(k_fcl_B)
        c_k_imp_1 = ctypes.c_double()
        c_k_imp_2 = ctypes.c_double()
        c_Omega_A = ctypes.c_double()
        c_Omega_B = ctypes.c_double()
        status = _two_objectives(c_k_fcl_A, c_k_fcl_B,
                                 ctypes.byref(c_k_imp_1),
                                 ctypes.byref(c_k_imp_2),
                                 ctypes.byref(c_Omega_A),
                                 ctypes.byref(c_Omega_B),
                                 ctypes.byref(self.pars),
                                 self._work)
        if status != 0:
            raise NumericException()

        return (c_k_imp_1.value, c_k_imp_2.value,
                c_Omega_A.value, c_Omega_B.value)

                            
