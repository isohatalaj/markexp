
import ctypes
import ctypes.util

max_banks = 2

# Try searching for the C library
bnetpath = ctypes.util.find_library("bnet")
if bnetpath:
    # Found, load it:
    libbnet = ctypes.CDLL(bnetpath)
else:
    # Not found, try loading Windows DLLs directly:
    try:
        # Try loading 32bit library
        libbnet = ctypes.CDLL("bnet32.dll")
    except OSError:
        # Didn't load. How about the 64bit version:
        libbnet = ctypes.CDLL("bnet64.dll")


_work_alloc = libbnet.bnet_work_alloc
_work_alloc.argtypes = []
_work_alloc.restype = ctypes.c_void_p

_work_free = libbnet.bnet_work_free
_work_free.argtypes = [ctypes.c_void_p]
_work_free.restype = None

_kzero = libbnet.bnet_kzero
_kzero.argtypes = [ctypes.c_double]
_kzero.restype = ctypes.c_double

_comp_kimp = libbnet.bnet_comp_kimp
_comp_kimp.argtypes = [ctypes.POINTER(ctypes.c_double),
                       ctypes.POINTER(ctypes.c_double),
                       ctypes.c_double, ctypes.c_double, ctypes.c_double]
_comp_kimp.restype = ctypes.c_double

_eps_star = libbnet.bnet_eps_star
_eps_star.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]
_eps_star.restype = ctypes.c_double

_comp_Omega = libbnet.bnet_comp_Omega
_comp_Omega.argtypes = [ctypes.c_double, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double,
                        ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double),
                        ctypes.POINTER(ctypes.c_double),
                        ctypes.c_void_p, ctypes.c_void_p]


class NumericException(Exception):
    pass

class Pars(ctypes.Structure):
    _fields_ = [
        ("LGD", ctypes.c_double),
        ("PD", ctypes.c_double),
        ("rho", ctypes.c_double),
        ("psi", ctypes.c_double),
        ("gamma", ctypes.c_double),
        ("q0", ctypes.c_double),
        ("q1", ctypes.c_double),
        ("phi", ctypes.c_double),
        ("phi2", ctypes.c_double),
        ("chi", ctypes.c_double),
        ("theta", ctypes.c_double)]


class Bnet:
    def __init__(self, *args, **kwargs):
        self._work = _work_alloc()
        self.pars = Pars(*args, **kwargs)

    def __del__(self):
        _work_free(self._work)

    def comp_kzero(self):
        return _kzero(self.pars.PD)

    def comp_epsstar(self, k, X):
        return _eps_star(k, X, self.pars.rho)

    def comp_kimp(self, epsfc, L, X):
        k0 = self.comp_kzero()
        c_epsfc = (ctypes.c_double * max_banks)(epsfc[0], epsfc[1])
        c_L = (ctypes.c_double * max_banks)(L[0], L[1])
        return _comp_kimp(ctypes.cast(c_epsfc, ctypes.POINTER(ctypes.c_double)),
                          ctypes.cast(c_L, ctypes.POINTER(ctypes.c_double)),
                          X, k0, self.pars.psi)


    def calibrate_phi(self):
        phical = ctypes.c_double()
        status = libbnet.bnet_calibrate_phi(self._work,
                                            ctypes.byref(self.pars),
                                            ctypes.byref(phical))
        if status != 0:
            raise NumericException()
        return phical.value

    def comp_Omegaplus(self, c0, X, kfc, kimp):
        Omega = ctypes.c_double()
        c1 = ctypes.c_double()
        Ec2 = ctypes.c_double()
        ipr = ctypes.c_double()
        c_c0 = ctypes.c_double(c0)
        c_X = ctypes.c_double(X)
        status = _comp_Omega(c0, X, kfc, kimp,
                             ctypes.byref(Omega),
                             ctypes.byref(c1),
                             ctypes.byref(Ec2),
                             ctypes.byref(ipr),
                             ctypes.byref(self.pars),
                             self._work)

        if status != 0:
            raise NumericException()

        return (Omega.value, c1.value, Ec2.value, ipr.value)




