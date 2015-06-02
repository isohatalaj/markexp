
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages 

import pybnet


# **********************************************************************
# Initialize solver with the given parameter values
bn = pybnet.Bnet(LGD = 0.2,
                 PD = 0.05,
                 rho = 0.3,
                 gamma = 0.2,
                 psi = 0.0,
                 q0 = 0.1,
                 q1 = 0.9,
                 chi = 0.0,
                 theta = 1.0)

# Calibrate phi, and set phi2 relative to phi
bn.pars.phi = bn.calibrate_phi()
bn.pars.phi2 = 2*bn.pars.phi

# Initial capitalisation ratio c, bank capital C, loan value L, and
# shock X.
c = (0.06, 0.04)
C = (4.5, 3.0)
L = (C[0]/c[0], C[1]/c[1])
X = -1.5

# Create vectorized versions for k impairment threshold and Omega_A
# and Omega_B.
K = np.vectorize(lambda ka, kb: 
                 bn.comp_kimp((bn.comp_epsstar(ka, X),
                               bn.comp_epsstar(kb, X)),
                              L, X))
WA = np.vectorize(lambda ka, kb, n: 
                  bn.comp_Omegaplus(C[0]/L[0], X, ka, K(ka, kb))[n])
WB = np.vectorize(lambda ka, kb, n: 
                  bn.comp_Omegaplus(C[1]/L[1], X, kb, K(ka, kb))[n])

# ---------------------------------------------------------------------
# NO PRICE INTERACTION
# Present parameters have psi = 0. Generate data.
kzero = bn.comp_kzero()    
kdata = np.linspace(kzero - 2.5, kzero, 40)

kbuse = kzero - 1.0;

G1data = WA(kdata, kbuse, 2)   
R1data = WA(kdata, kbuse, 3)

# ---------------------------------------------------------------------
# STRONG PRICE INTERACTION
# Modify parameters, re-calibrate (unnecessary?), and generate data.
bn.pars.psi = 7.0
bn.pars.phi = bn.calibrate_phi()
bn.pars.phi2 = 2*bn.pars.phi

G2data = WA(kdata, kbuse, 2)
R2data = WA(kdata, kbuse, 3)

# ---------------------------------------------------------------------
# PLOT BOTH TOGETHER

mpl.rc('text', usetex = True)
mpl.rc('font', size=16)

gaincol = '#0072bd'
riskcol = '#d95319'

fig = plt.figure(figsize=(5, 4))

plt.plot(kdata, G2data, '-', label=r'$\mbox{E}[C^2], \psi = 7$', color=gaincol)
plt.plot(kdata, G1data, '--', label=r'$\mbox{E}[C^2], \psi = 0$', color=gaincol)
plt.plot(kdata, R2data, '-', label=r'$\mbox{IPR}_{80}[C^2 - C^1], \psi = 7$', color=riskcol)
plt.plot(kdata, R1data, '--', label=r'$\mbox{IPR}_{80}[C^2 - C^1], \psi = 0$', color=riskcol)

plt.xlabel(r'$\longleftarrow \mbox{Refinance} \qquad k^{\dag\dag}_A \qquad \mbox{Foreclose} \longrightarrow$')
plt.ylabel(r'$\mbox{E}[C^2], \quad \mbox{IPR}_{80}[C^2 - C^1]$')

plt.legend(bbox_to_anchor=(0,0), loc=3, frameon=False, prop={'size': 12})

plt.subplots_adjust(left=0.175, bottom=0.175, 
                    right=0.95, top=0.95)

# Set here output to None to get output to screen.
# output = None
output = 'markexp-fig1.pdf'

if output:
    pp = PdfPages(output)
    pp.savefig()
    pp.close()
else:
    plt.show()

plt.close()

