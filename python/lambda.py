
import math
import numpy as np
import basicplot3d as p3d

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages 


def reval1(eps1, epsimp1, epsfc1, phi):
    if eps1 < epsfc1:
        return 0.0
        
    if eps1 > epsimp1:
        return 1.0

    return math.exp(-phi*(epsimp1 - eps1))

def reval2(eps1, eps2, epsimp1, epsfc1, epsimp2, phi, varphi):
    theta = eps2 - epsimp2 - varphi*max(epsimp1-eps1, 0.0)

    if eps1 < epsfc1:
        return 0.0

    if theta > 0.0:
        return 1.0

    return math.exp(phi*theta)

def plot_reval1(epsimp1, epsfc1, phi):
    mpl.rc('text', usetex = True)
    mpl.rc('font', size=24)

    e1data = np.linspace(-3.0, 1.0, 200)

    lambdafun = np.vectorize(lambda e1:
                             reval1(e1, epsimp1, epsfc1, phi))

    ldata = lambdafun(e1data)
    plt.plot(e1data, ldata)

    plt.tight_layout(pad=2)

    plt.axis([-3.0, 1.0, 0.0, 1.33])
    plt.xlabel(r'$\epsilon_1^i$')
    plt.ylabel(r'$\lambda_1^i$')

    plt

    # output = None
    output = 'lambda1.pdf'

    
    if output:
        pp = PdfPages(output)
        pp.savefig()
        pp.close()
    else:
        plt.show()
        
    plt.close()


def plot_reval2(epsimp1, epsfc1, epsimp2, phi, varphi):
    e1data = np.linspace(-3.0, 1.0, 30)
    e2data = np.linspace(-3.0, 1.0, 30)
    e1data, e2data = np.meshgrid(e1data, e2data)

    lambdafun = np.vectorize(lambda e1, e2:
                             reval2(e1, e2, epsimp1, epsfc1, epsimp2, 
                                   phi, varphi))

    ldata = lambdafun(e1data, e2data)

    p3d.plot3d(e1data, e2data, ldata,
               xlabel=r'$\epsilon_1^i$',
               ylabel=r'$\epsilon_2^i$',
               zlabel=r'$\lambda^i_2$'
               , export="lambda2.pdf"
    )


plot_reval2(-1.0, -2.0, -1.5, 1.0, 1.0)
plot_reval1(-1.0, -2.0, 1.0)

               
    



