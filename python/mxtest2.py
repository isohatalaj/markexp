

import numpy as np
import mayavi.mlab as my

from traits.api import HasTraits, Range, Instance, \
        on_trait_change
from traitsui.api import View, Item, Group, HSplit, VSplit

from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import MayaviScene, SceneEditor, \
                MlabSceneModel


# Load the library containing an interface to our C solver
import pymexp


# **********************************************************************
# INITIALIZATION

class MEModel(HasTraits):
    _phi = Range(0.01, 2.0, 0.5)
    _varphi = Range(0.01, 2.0, 1.0)
    _k0 = Range(-4.0, 0.0, -1.5)
    _L0_A = Range(0.0, 100.0, 20.0)
    _L0_B = Range(0.0, 100.0, 40.0)
    _rho = Range(0.01, 0.99, 0.3)
    _xi = Range(0.01, 0.99, 0.3)
    _mu = Range(0.0, 1.0, 0.0)
    _zeta = Range(0.0, 1.0, 1.0)
    _psi = Range(0.0, 10.0, 1.0)
    _chi = Range(0.0, 1.0, 1.0)
    _gamma = Range(0.0, 1.0, 0.5)
    _q0 = 0.1
    _q1 = 0.9
    _X1 = Range(-3.0, 3.0, -1.0)
    
    
    # Initial capitalisation ratio c, bank capital C, loan value L, and
    # shock X.
    # c = (0.06, 0.06)
    # C = (4.5, 3.0)
    # L = (C[0]/c[0], C[1]/c[1])
    # X = -1.5

    first = True

    scenea = Instance(MlabSceneModel, ())
    sceneb = Instance(MlabSceneModel, ())
    scenec = Instance(MlabSceneModel, ())

    surfa = Instance(PipelineBase)
    surfb = Instance(PipelineBase)
    surfc = Instance(PipelineBase)

    def __init__(self):
        super(HasTraits, self).__init__()

        # Initialize solver with the given parameter values
        self.me = pymexp.Mexp(phi1 = self._phi,
                              phi2 = self._phi,
                              varphi = self._varphi,
                              k0_1 = self._k0,
                              k0_2 = self._k0,
                              L0_A = self._L0_A,
                              L0_B = self._L0_B,
                              rho = self._rho,
                              xi = self._xi,
                              mu = self._mu,
                              zeta = self._zeta,
                              psi = self._psi,
                              chi = self._chi,
                              gamma = self._gamma,
                              q0 = self._q0,
                              q1 = self._q1,
                              X1 = self._X1)


    @on_trait_change('scenec.activated')
    def populate_scenes(self):
        self.update_plot()

    @on_trait_change('_phi, _varphi, _k0, _L0_A, _L0_B, _rho, _xi, _mu, _zeta, _psi, _chi, _gamma, _X1')
    def update_plot(self):

        # Disable rendering for the duration of the update
        self.scenea.disable_render = True
        self.sceneb.disable_render = True
        self.scenec.disable_render = True

        # Read parameters from the UI
        self.me.pars.phi1 = self._phi
        self.me.pars.phi2 = self._phi
        self.me.pars.varphi = self._varphi
        self.me.pars.k0_1 = self._k0
        self.me.pars.k0_2 = self._k0
        self.me.pars.L0_A = self._L0_A
        self.me.pars.L0_B = self._L0_B
        self.me.pars.rho = self._rho
        self.me.pars.xi = self._xi
        self.me.pars.mu = self._mu
        self.me.pars.zeta = self._zeta
        self.me.pars.psi = self._psi
        self.me.pars.chi = self._chi
        self.me.pars.gamma = self._gamma
        self.me.pars.q0 = self._q0
        self.me.pars.q1 = self._q1
        self.me.pars.X1 = self._X1

        xmin = self.me.pars.k0_1 - 3.0
        xmax = self.me.pars.k0_1 
        ymin = self.me.pars.k0_1 - 3.0
        ymax = self.me.pars.k0_1 
        d = 8*1j
        
        xdata, ydata = np.mgrid[xmin:xmax:d, ymin:ymax:d]

        self.F = np.vectorize(lambda k_fcl_A, k_fcl_B, n: 
                              self.me.objectives(k_fcl_A, k_fcl_B)[n])

        z0data = self.F(xdata, ydata, 2)
        z1data = self.F(xdata, ydata, 3)
        z2data = z1data + z0data

        z0 = z0data.min()
        z1 = z0data.max()
        z0data = (z0data - z0) / (z1 - z0)

        z0 = z1data.min()
        z1 = z1data.max()
        z1data = (z1data - z0) / (z1 - z0)

        z0 = z2data.min()
        z1 = z2data.max()
        z2data = (z2data - z0) / (z1 - z0)


        if self.first:
            self.first = False
            self.surfa = self.scenea.mlab.surf(xdata, ydata, z0data, figure=self.scenea.mayavi_scene)#, warp_scale='auto')
            self.surfb = self.sceneb.mlab.surf(xdata, ydata, z1data, figure=self.sceneb.mayavi_scene)#, warp_scale='auto')
            self.surfc = self.scenec.mlab.surf(xdata, ydata, z2data, figure=self.scenec.mayavi_scene)#, warp_scale='auto')

            my.axes(self.surfa, nb_labels=5, xlabel="k_fcl_A", ylabel="k_fcl_B", zlabel="Omega_A", figure=self.scenea.mayavi_scene)
            my.axes(self.surfb, nb_labels=5, xlabel="k_fcl_A", ylabel="k_fcl_B", zlabel="Omega_B", figure=self.sceneb.mayavi_scene)
            my.axes(self.surfc, nb_labels=5, xlabel="k_fcl_A", ylabel="k_fcl_B", zlabel="Omega_A+B", figure=self.scenec.mayavi_scene)

            my.sync_camera(self.scenea.mayavi_scene, self.sceneb.mayavi_scene)
            my.sync_camera(self.scenea.mayavi_scene, self.scenec.mayavi_scene)

        else:
            self.surfa.mlab_source.set(x=xdata, y=ydata, scalars=z0data)
            self.surfb.mlab_source.set(x=xdata, y=ydata, scalars=z1data)
            self.surfc.mlab_source.set(x=xdata, y=ydata, scalars=z2data)

        # Resume rendering
        self.scenea.disable_render = False
        self.sceneb.disable_render = False
        self.scenec.disable_render = False

        # Update done.


    # The layout of the dialog created
    view = View(VSplit(HSplit(Item('scenea', editor=SceneEditor(scene_class=MayaviScene),
                                   height=250, width=300, show_label=False),
                              Item('sceneb', editor=SceneEditor(scene_class=MayaviScene),
                                   height=250, width=300, show_label=False),
                              Item('scenec', editor=SceneEditor(scene_class=MayaviScene),
                                   height=250, width=300, show_label=False),
                              ),
                       Group('_', '_phi', '_varphi', '_k0',
                             '_L0_A', '_L0_B', '_rho', '_xi',
                             '_mu', '_zeta', '_psi', '_chi',
                             '_gamma', '_X1'),
                       ),
                resizable=True,
                )

memod = MEModel()
memod.configure_traits()
