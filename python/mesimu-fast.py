

import numpy as np
import mayavi.mlab as my

# Load the library containing an interface to our C solver
import pybnet


# **********************************************************************
# INITIALIZATION

# Initialize solver with the given parameter values
bn = pybnet.Bnet(LGD = 0.2,
                 PD = 0.05,
                 rho = 0.3,
                 gamma = 0.5,
                 psi = 5.0,
                 q0 = 0.1,
                 q1 = 0.9,
                 chi = 1.0)

# Calibrate phi, and phi2 relative to phi
phi_calib = bn.calibrate_phi()

bn.pars.phi = phi_calib
bn.pars.phi2 = 2*bn.pars.phi

# Initial capitalisation ratio c, bank capital C, loan value L, and
# shock X.
c = (0.06, 0.045)
C = (4.5, 3.0)
L = (C[0]/c[0], C[1]/c[1])
X = -1.5


# **********************************************************************

# For conveniance, we will make vectorized versions versions of of
# the kimp and OmegaA, OmegaB functions. All this does is it allows
# us to apply these functions directly to arrays of points, rather
# than having to write explicit loops.

K = np.vectorize(lambda ka, kb: 
                 bn.comp_kimp((bn.comp_epsstar(ka, X),
                               bn.comp_epsstar(kb, X)),
                              L, X))
WA = np.vectorize(lambda ka, kb, n: 
                  bn.comp_Omegaplus(C[0]/L[0], X, ka, K(ka, kb))[n])
WB = np.vectorize(lambda ka, kb, n: 
                  bn.comp_Omegaplus(C[1]/L[1], X, kb, K(ka, kb))[n])


# **********************************************************************
# COMPUTE DATA

# Let us choose foreclosure thresholds (kfcA, kfcB) relative to
# baseline k, kzero;
kzero = bn.comp_kzero()

# kadata = np.ogrid[kzero - 2.5:kzero:20j]
# kbdata = np.ogrid[kzero - 2.5:kzero:20j]
kadata, kbdata = np.mgrid[kzero - 2.5:kzero:20j, kzero - 2.5:kzero:20j]

# kadata = np.linspace(kzero - 2.5, kzero, 20)
# kbdata = np.linspace(kzero - 2.5, kzero, 20)
# kadata, kbdata = np.meshgrid(kadata, kbdata)

# WAdata = WA(kadata, kbdata, 0)
# WBdata = WB(kadata, kbdata, 0)
# WSdata = WAdata + WBdata

# WA(kadata, kbdata, 1) would give bank A C^1 vs. foreclosure k
# WA(kadata, kbdata, 2) would give bank A E[C^2] vs. foreclosure k
# WA(kadata, kbdata, 3) would give bank A IPR[C^2-C^1] vs. foreclosure k
# WB gives the above for bank B

# **********************************************************************
# PLOT DATA

figa = my.figure()
figb = my.figure()
figc = my.figure()

def draw():
    surf = my.surf(kadata, kbdata, (lambda ka, kb: WA(ka, kb, 0)),
                   warp_scale='auto',
                   figure=figa)
    my.axes(surf, nb_labels=5)
    
    surf = my.surf(kadata, kbdata, (lambda ka, kb: WB(ka, kb, 0)),
                   warp_scale='auto',
                   figure=figb)
    my.axes(surf, nb_labels=5)
    
    surf = my.surf(kadata, kbdata, (lambda ka, kb: WA(ka, kb, 0) + WB(ka, kb, 0)),
                   warp_scale='auto',
                   figure=figc)
    my.axes(surf, nb_labels=5)
    
    my.sync_camera(figa, figb)
    my.sync_camera(figa, figc)

draw()

my.show()


# **********************************************************************
# CLEANUP
# Destroy the solver object.
# del bn
