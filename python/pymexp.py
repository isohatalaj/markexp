
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
_two_objectives.restype = ctypes.c_int

_find_k_fcl_global_ubound = libmexp.mexp_find_k_fcl_global_ubound
_find_k_fcl_global_ubound.argtypes = [ctypes.POINTER(ctypes.c_double),
                               ctypes.c_void_p, ctypes.c_void_p]
_find_k_fcl_global_ubound.restype = ctypes.c_int

_k_fcls_from_ps = libmexp.mexp_k_fcls_from_ps
_k_fcls_from_ps.argtypes = [ctypes.c_double, ctypes.c_double,
                            ctypes.POINTER(ctypes.c_double),
                            ctypes.POINTER(ctypes.c_double),
                            ctypes.c_void_p, ctypes.c_void_p]
_k_fcls_from_ps.restype = ctypes.c_int

_eps_star = libmexp.mexp_eps_star
_eps_star.argtypes = [ctypes.c_double, ctypes.c_void_p]
_eps_star.restype = ctypes.c_double

_N_star = libmexp.mexp_N_star
_N_star.argtypes = [ctypes.c_double, ctypes.c_void_p]
_N_star.restype = ctypes.c_double

_k_of_p_tilde = libmexp.mexp_k_of_p_tilde
_k_of_p_tilde.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_void_p]
_k_of_p_tilde.restype = ctypes.c_double

_cap2_pdf_cdf = libmexp.mexp_cap2_pdf_cdf
_cap2_pdf_cdf.argtypes = [ctypes.c_double, ctypes.c_double, 
                          ctypes.c_int, ctypes.c_int,
                          ctypes.c_double, 
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.POINTER(ctypes.c_double),
                          ctypes.c_void_p, ctypes.c_void_p]
_cap2_pdf_cdf.restype = ctypes.c_int

class NumericException(Exception):
    pass
                 
class Pars(ctypes.Structure):
    _fields_ = [
        ("phi1", ctypes.c_double),
        ("phi2", ctypes.c_double),
        ("varphi", ctypes.c_double),
        ("k0_1", ctypes.c_double),
        ("k0_2", ctypes.c_double),
        ("c0_A", ctypes.c_double),
        ("c0_B", ctypes.c_double),
        ("L0_A", ctypes.c_double),
        ("L0_B", ctypes.c_double),
        ("rho", ctypes.c_double),
        ("xi", ctypes.c_double),
        ("mu", ctypes.c_double),
        ("zeta", ctypes.c_double),
        ("psi", ctypes.c_double),
        ("chi", ctypes.c_double),
        ("gamma0", ctypes.c_double),
        ("gamma1", ctypes.c_double),
        ("q0", ctypes.c_double),
        ("q1", ctypes.c_double),
        ("gamma2", ctypes.c_double),
        ("c_bar", ctypes.c_double),
        ("X1", ctypes.c_double),
        ("dummy", ctypes.c_double)
    ]

class Mexp:
    def __init__(self, *args, **kwargs):
        self._work = _work_alloc()
        self.pars = Pars(*args, **kwargs)

    def __del__(self):
        _work_free(self._work)

    def set_baseline(self):
        phi = 0.2
        k0 = -1.5
        c0 = 0.07
        self.pars.phi1 = phi
        self.pars.phi2 = phi
        self.pars.varphi = 1.0
        self.pars.k0_1 = k0
        self.pars.k0_2 = k0
        self.pars.c0_A = c0
        self.pars.c0_B = c0
        self.pars.L0_A = 50
        self.pars.L0_B = 50
        self.pars.rho = 0.4
        self.pars.xi = 0.6
        self.pars.mu = 0.5
        self.pars.zeta = 0.866
        self.pars.psi = 0.7
        self.pars.chi = 0.0
        self.pars.gamma0 = 0.25
        self.pars.gamma1 = 0.0
        self.pars.gamma2 = 1.0
        self.pars.c_bar = 0.04
        self.pars.q0 = 0.1
        self.pars.q1 = 0.9
        self.pars.X1 = -2.0

    # def map_coefs(self, eps_imp_1, eps_fcl_1, eps_imp_2):
    #     c_eps_imp_1 = ctypes.c_double(eps_imp_1)
    #     c_eps_fcl_1 = ctypes.c_double(eps_fcl_1)
    #     c_eps_imp_2 = ctypes.c_double(eps_imp_2)
    #     c_Lambda_1 = ctypes.c_double()
    #     c_Kappa_1 = ctypes.c_double()
    #     c_Iota_1 = ctypes.c_double()
    #     c_Lambda_2 = ctypes.c_double()
    #     status = _map_coefs(c_eps_imp_1, c_eps_fcl_1, c_eps_imp_2,
    #                         ctypes.byref(c_Iota_1),
    #                         ctypes.byref(c_Kappa_1),
    #                         ctypes.byref(c_Lambda_1),
    #                         ctypes.byref(c_Lambda_2),
    #                         ctypes.byref(self.pars),
    #                         self._work)
    #     if status != 0:
    #         raise NumericException()

    #     return (c_Iota_1.value, c_Kappa_1.value, c_Lambda_1.value, c_Lambda_2.value)

    # def comp_P1(self, eps_imp_1, eps_fcl_A_1, eps_fcl_B_1):
    #     (IoA, KpA, Lm1A, Lm2A) = self.map_coefs(eps_imp_1, eps_fcl_A_1, 0.0)
    #     (IoB, KpB, Lm1B, Lm2B) = self.map_coefs(eps_imp_1, eps_fcl_B_1, 0.0)
        
    #     return ((self.pars.c0_A + Kp)*self.pars.L0_A, (self.pars.c0_A + Kp)/Lm1,
    #             (self.pars.c0_B + Kp)*self.pars.L0_B, (self.pars.c0_B + Kp)/Lm1)

    # def comp_cr2(self, eps_imp_1, eps_fcl_1, eps_imp_2):
    #     (Io, Kp, Lm1, Lm2) = self.map_coefs(eps_imp_1, eps_fcl_1, eps_imp_2)
    #     return (self.pars.c0_A + Kp, (self.pars.c0_A + Kp)/Lm1,
    #             self.pars.c0_A + Io - 1 + Lm2, (self.pars.c0_A + Io - 1)/Lm2 + 1.0)

    def find_k_fcl_global_ubound(self):
        c_k_fcl_ubound = ctypes.c_double()
        status = _find_k_fcl_global_ubound(ctypes.byref(c_k_fcl_ubound),
                                           ctypes.byref(self.pars),
                                           self._work)
        if status != 0:
            raise NumericException()
        
        return c_k_fcl_ubound.value

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

    # bank = 0 for A, 1 for B
    # fun = 0 for absolute C, 1 for cap ratio
    def cap2_pdf_cdf(self, bank, fun, x, k_fcl_A, k_fcl_B):
        c_k_fcl_A = ctypes.c_double(k_fcl_A)
        c_k_fcl_B = ctypes.c_double(k_fcl_B)
        c_x = ctypes.c_double(x)
        c_pdf = ctypes.c_double()
        c_cdf = ctypes.c_double()
        c_fun = ctypes.c_int(fun)
        c_bank = ctypes.c_int(bank)
        status = _cap2_pdf_cdf(c_k_fcl_A, c_k_fcl_B,
                               c_bank, c_fun,
                               c_x, 
                               ctypes.byref(c_pdf),
                               ctypes.byref(c_cdf),
                               ctypes.byref(self.pars),
                               self._work)
        if status != 0:
            raise NumericException()

        return (c_pdf.value, c_cdf.value)
        
    def k_fcls_from_ps(self, p_A, p_B):
        c_p_A = ctypes.c_double(p_A)
        c_p_B = ctypes.c_double(p_B)
        c_k_fcl_A = ctypes.c_double()
        c_k_fcl_B = ctypes.c_double()

        status = _k_fcls_from_ps(c_p_A, c_p_B,
                                 ctypes.byref(c_k_fcl_A),
                                 ctypes.byref(c_k_fcl_B),
                                 ctypes.byref(self.pars),
                                 self._work)
        
        if status != 0:
            raise NumericException()

        return (c_k_fcl_A.value, c_k_fcl_B.value)
                                 
    def objectives_p(self, p_A, p_B):
        k_fcl_A, k_fcl_B = self.k_fcls_from_ps(p_A, p_B)
        return self.objectives(k_fcl_A, k_fcl_B)


    # def c2_A_pdf_cdf(self, c2, k_fcl_A, k_fcl_B): 
    #     c_k_fcl_A = ctypes.c_double(k_fcl_A)
    #     c_k_fcl_B = ctypes.c_double(k_fcl_B)
    #     c_c2 = ctypes.c_double(c2)
    #     c_pdf = ctypes.c_double()
    #     c_cdf = ctypes.c_double()
    #     status = _c2_A_pdf_cdf(c_k_fcl_A, c_k_fcl_B,
    #                            c_c2, 
    #                            ctypes.byref(c_pdf),
    #                            ctypes.byref(c_cdf),
    #                            ctypes.byref(self.pars),
    #                            self._work)
    #     if status != 0:
    #         raise NumericException()

    #     return (c_pdf.value, c_cdf.value)


    def eps_star(self, k):
        c_k = ctypes.c_double(k)
        return _eps_star(k, ctypes.byref(self.pars))

    def N_star(self, k):
        c_k = ctypes.c_double(k)
        return _N_star(k, ctypes.byref(self.pars))

    def k_of_p_tilde(self, p_tilde, k_bound):
        return _k_of_p_tilde(ctypes.c_double(p_tilde),
                             ctypes.c_double(k_bound),
                             ctypes.byref(self.pars))
