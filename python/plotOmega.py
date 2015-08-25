
import numpy as np

import pymexp
import basicplot3d as p3d

me = pymexp.Mexp()
me.set_baseline()

k_max = me.pars.k0_1
d = 32j

xdata, ydata = np.mgrid[0.01:0.99:d, 0.01:0.99:d]
Omega = np.vectorize(lambda x, y, n: 
                     me.objectives(me.k_of_p_tilde(x, k_max),
                                   me.k_of_p_tilde(y, k_max))[n])

OmegaAdata = Omega(xdata, ydata, 2)
OmegaBdata = Omega(xdata, ydata, 3)
OmegaCdata = OmegaAdata + OmegaBdata

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
           export='Omega-A.pdf',
           **kwargs)

p3d.plot3d(xdata, ydata, OmegaBdata,
           zlabel=r'$\Omega_B$',
           plotlabel=r'$(b)$',
           export='Omega-B.pdf',
           **kwargs)

p3d.plot3d(xdata, ydata, OmegaCdata,
           zlabel=r'$\Omega_A + \Omega_B$',
           plotlabel=r'$(c)$',
           export='Omega-C.pdf',
           **kwargs)

