
import numpy as np

import pymexp
import basicplot3d as p3d

# ************************************************************************
# Setup solver

me = pymexp.Mexp()
me.set_baseline()


# ************************************************************************
# Generate data

dx = 128j
dy = 24j

c2_lo = 0.031
c2_hi = 0.069

C2_lo = 0.8*me.pars.L0_A*c2_lo
C2_hi = 1.0*me.pars.L0_A*c2_hi

x1data, y1data = np.mgrid[c2_lo:c2_hi:dx, 0.01:0.99:dy]
x2data, y2data = np.mgrid[C2_lo:C2_hi:dx, 0.01:0.99:dy]

p_B = 0.5

fF = np.vectorize(lambda x, y, b, f, n:
                      me.cap2_pdf_cdf_p(b, f, x, y, p_B)[n])

c2Apdf = fF(x1data, y1data, 0, 1, 0)
C2Apdf = fF(x2data, y2data, 0, 0, 0)

# ************************************************************************
# Additional info to screen

print "P(c2 < cmin) = ", fF(me.pars.c_bar, 0.01, 0, 1, 1), "(pA = 0.01) ...",
print fF(me.pars.c_bar, 0.99, 0, 1, 1), "(pA = 0.99)"


# ************************************************************************
# Setup plot parameters

# 3D view angles, rotation around z axis (azimuth) and elevation above
# xy-plane:
azim = -60
elev = 55

# Trim output so that the graph is neatly contained in the figure
# area. Vector of left, bottom, right, and top trim values. Negative
# shrinks, positive enlarges in the corresponding direction.
trim = [-0.1, 0.1, 0.1, 0.0]

# Package the plot options
kwargs = {"elev": elev, "azim": azim, "trim": trim,
          "zform": "5.0f", "zlabelshift": 3, 
          "labeldist": 3}


# ************************************************************************
# Do plots

p3d.plot3d(x1data, y1data, c2Apdf,
           xlabel=r'$C_2/L_2$',
           ylabel=r'$p_A$',
           zlabel=r'$C_2/L_2\; \mathrm{pdf.}$',
           plotlabel="(a)",
           export='caprat2-pdf-A.pdf',
           **kwargs)

p3d.plot3d(x2data, y2data, C2Apdf,
           xlabel=r'$C_2$',
           ylabel=r'$p_A$',
           zlabel=r'$C_2\; \mathrm{pdf.}$',
           plotlabel="(b)",
           export='caplev2-pdf-A.pdf',
           **kwargs)

