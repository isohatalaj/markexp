
import numpy as np

import pybnet
import pybnplots as plot

def plot_Omega(K, WA, WB, fileprefix):
    """
    Plot Omega_A, Omega_B, and Omega_A + Omega_B over the k
    foreclosure thresholds, given the price adjusted impairment
    threshold function K. Parameters WA and WB are vectorized Omega_A
    and Omega_B functions.
    """

    kzero = bn.comp_kzero()    

    kadata = np.linspace(kzero - 2.5, kzero, 20)
    kbdata = np.linspace(kzero - 2.5, kzero, 20)
    kadata, kbdata = np.meshgrid(kadata, kbdata)
    
    WAdata = WA(kadata, kbdata, 0)
    WBdata = WB(kadata, kbdata, 0)
    WSdata = WAdata + WBdata

    plot.plot3d(kadata, kbdata, 1000*WAdata,
                xlabel=r'$k^{\dag\dag}_A$',
                ylabel=r'$k^{\dag\dag}_B$',
                zlabel=r'$\Omega_A \quad [1/1000]$',
                plotlabel=r'$(a)$',
                export=fileprefix+'Omega-A.pdf')

    plot.plot3d(kadata, kbdata, 1000*WBdata,
                xlabel=r'$k^{\dag\dag}_A$',
                ylabel=r'$k^{\dag\dag}_B$',
                zlabel=r'$\Omega_B \quad [1/1000]$',
                plotlabel=r'$(b)$',
                export=fileprefix+'Omega-B.pdf')

    plot.plot3d(kadata, kbdata, 1000*WSdata,
                xlabel=r'$k^{\dag\dag}_A$',
                ylabel=r'$k^{\dag\dag}_B$',
                zlabel=r'$\Omega_A + \Omega_B \quad [1/1000]$',
                plotlabel=r'$(c)$',
                export=fileprefix+'Omega-A+B.pdf')


# **********************************************************************
# FIGURE 2

# Initialize solver with the given parameter values
bn = pybnet.Bnet(LGD = 0.2,
                 PD = 0.05,
                 rho = 0.3,
                 gamma = 0.2,
                 psi = 1.0,
                 q0 = 0.1,
                 q1 = 0.9,
                 chi = 1.0)

# Calibrate phi, and phi2 relative to phi
phi_calib = bn.calibrate_phi()

bn.pars.phi = phi_calib
bn.pars.phi2 = 2*bn.pars.phi

# Initial capitalisation ratio c, bank capital C, loan value L, and
# shock X.
c = (0.06, 0.04)
C = (4.5, 3.0)
L = (C[0]/c[0], C[1]/c[1])
X = -1.5


K = np.vectorize(lambda ka, kb: 
                 bn.comp_kimp((bn.comp_epsstar(ka, X),
                               bn.comp_epsstar(kb, X)),
                              L, X))
WA = np.vectorize(lambda ka, kb, n: 
                  bn.comp_Omegaplus(C[0]/L[0], X, ka, K(ka, kb))[n])
WB = np.vectorize(lambda ka, kb, n: 
                  bn.comp_Omegaplus(C[1]/L[1], X, kb, K(ka, kb))[n])


plot_Omega(K, WA, WB, "markexp-fig2-")


# **********************************************************************
# FIGURE 2

# Strong risk weight, strong price interaction
bn.pars.gamma = 0.5
bn.pars.psi = 5.0

bn.pars.phi = bn.calibrate_phi()
bn.pars.phi2 = 2*bn.pars.phi

# The functions K, WA, WB remain bound to the bn object, so we need
# not recompute them.

plot_Omega(K, WA, WB, "markexp-fig3-")
