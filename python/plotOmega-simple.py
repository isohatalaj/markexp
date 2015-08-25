
import numpy as np

import pymexp
import basicplot3d as p3d

# ************************************************************************
# Setup solver
me = pymexp.Mexp()
me.set_baseline()

# ************************************************************************
k_max = me.pars.k0_1
d = 32j

xdata, ydata = np.mgrid[0.01:0.99:d, 0.01:0.99:d]
Omega = np.vectorize(lambda x, y, n: 
                     me.objectives(me.k_of_p_tilde(x, k_max),
                                   me.k_of_p_tilde(y, k_max))[n])

me.pars.gamma0 = 1.0
me.pars.gamma1 = 0.0
me.pars.gamma2 = 0.0
OmegaAdata = Omega(xdata, ydata, 2)

me.pars.gamma0 = 0.0
me.pars.gamma1 = 1.0
me.pars.gamma2 = 0.0
OmegaBdata = Omega(xdata, ydata, 2)

me.pars.gamma0 = 0.0
me.pars.gamma1 = 0.0
me.pars.gamma2 = 1.0
OmegaCdata = Omega(xdata, ydata, 2)

# 3D view angles, rotation around z axis (azimuth) and elevation above
# xy-plane:
elev = 45
azim = -70

# Trim output so that the graph is neatly contained in the figure
# area. Vector of left, bottom, right, and top trim values. Negative
# shrinks, positive enlarges in the corresponding direction.
trim = [-0.1, 0.1, 0.1, 0.0]

# Package the plot options
kwargs = {"xlabel": r'$p_A$',
          "ylabel": r'$p_B$',
          "elev": elev, "azim": azim, "trim": trim,
          "zform": "5.2f", "zlabelshift": 7, 
          "labeldist": 3}

p3d.plot3d(xdata, ydata, OmegaAdata,
           zlabel=r'$\Omega_A$',
           plotlabel=r'$(a)$',
           export='Omega-simple-A.pdf',
           **kwargs)

p3d.plot3d(xdata, ydata, OmegaBdata,
           zlabel=r'$\Omega_A$',
           plotlabel=r'$(b)$',
           export='Omega-simple-B.pdf',
           **kwargs)

p3d.plot3d(xdata, ydata, OmegaCdata,
           zlabel=r'$\Omega_A$',
           plotlabel=r'$(c)$',
           export='Omega-simple-C.pdf',
           **kwargs)
