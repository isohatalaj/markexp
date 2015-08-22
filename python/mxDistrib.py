

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
    _phi = Range(0.01, 1.0, 0.5)
    _varphi = Range(0.01, 2.0, 1.0)
    _k0 = Range(-4.0, 0.0, -1.5)
    _c0 = Range(0.0, 0.25, 0.12)
    _L0_A = Range(0.0, 100.0, 1.0)
    _L0_B = Range(0.0, 100.0, 1.0)
    _rho = Range(0.01, 0.99, 0.3)
    _xi = Range(0.01, 0.99, 0.5)
    _mu = Range(0.0, 1.0, 0.45)
    _zeta = Range(0.0, 10.0, 1.0)
    _psi = Range(0.0, 10.0, 5.0)
    _chi = Range(0.0, 1.0, 0.0)
    _gamma0 = Range(0.0, 1.0, 1.0)
    _gamma1 = Range(0.0, 1.0, 0.5)
    _gamma2 = Range(0.0, 1.0, 0.0)
    _q0 = 0.1
    _q1 = 0.9
    _X1 = Range(-3.0, 3.0, -2.0)
    _c_bar = Range(-1.0, 0.1, 0.04)
    _ptilde_B = Range(0.0, 1.0, 0.5)
    _mesh_d = Range(4, 64, 12, )
    
    
    # Initial capitalisation ratio c, bank capital C, loan value L, and
    # shock X.
    # c = (0.06, 0.06)
    # C = (4.5, 3.0)
    # L = (C[0]/c[0], C[1]/c[1])
    # X = -1.5

    first = True

    scenea = Instance(MlabSceneModel, ())
    sceneb = Instance(MlabSceneModel, ())

    surfa = Instance(PipelineBase)
    surfb = Instance(PipelineBase)


    def __init__(self):
        super(HasTraits, self).__init__()

        # Initialize solver with the given parameter values
        self.me = pymexp.Mexp(phi1 = self._phi,
                              phi2 = self._phi,
                              varphi = self._varphi,
                              k0_1 = self._k0,
                              k0_2 = self._k0,
                              c0_A = self._c0,
                              c0_B = self._c0,
                              L0_A = self._L0_A,
                              L0_B = self._L0_B,
                              rho = self._rho,
                              xi = self._xi,
                              mu = self._mu,
                              zeta = self._zeta,
                              psi = self._psi,
                              chi = self._chi,
                              gamma0 = self._gamma0,
                              gamma1 = self._gamma1,
                              gamma2 = self._gamma2,
                              q0 = self._q0,
                              q1 = self._q1,
                              X1 = self._X1)


    @on_trait_change('sceneb.activated')
    def populate_scenes(self):
        self.update_plot()

    @on_trait_change('_phi, _varphi, _k0, _c0, _L0_A, _L0_B, _rho, _xi, _mu, _zeta, _psi, _chi, _gamma0, _gamma1, _gamma2, _X1, _c0, _c_bar, _mesh_d, _ptilde_B')
    def update_plot(self):
        
        print "-----------------------------------------------------------"
        # Disable rendering for the duration of the update
        self.scenea.disable_render = True
        self.sceneb.disable_render = True

        # Read parameters from the UI
        self.me.pars.phi1 = self._phi
        self.me.pars.phi2 = self._phi
        self.me.pars.varphi = self._varphi
        self.me.pars.k0_1 = self._k0
        self.me.pars.k0_2 = self._k0
        self.me.pars.c0_A = self._c0
        self.me.pars.c0_B = self._c0
        self.me.pars.L0_A = self._L0_A
        self.me.pars.L0_B = self._L0_B
        self.me.pars.rho = self._rho
        self.me.pars.xi = self._xi
        self.me.pars.mu = self._mu
        self.me.pars.zeta = self._zeta
        self.me.pars.psi = self._psi
        self.me.pars.chi = self._chi
        self.me.pars.gamma0 = self._gamma0
        self.me.pars.gamma1 = self._gamma1
        self.me.pars.gamma2 = self._gamma2
        self.me.pars.q0 = self._q0
        self.me.pars.q1 = self._q1
        self.me.pars.c_bar = self._c_bar;
        self.me.pars.X1 = self._X1

        k_fcl_ubound = self.me.find_k_fcl_ubound()
        print "k_fcl_ubound = ", k_fcl_ubound,
        print " => max fcl-% = ", self.me.N_star(k_fcl_ubound), 
        print "(baseline ", self.me.N_star(self.me.pars.k0_1), ")"


        ## k_max = k_fcl_ubound
        k_max = self.me.pars.k0_1

        d = self._mesh_d*1j
        
        # xmin = k_max - 3.0
        # xmax = k_max
        # ymin = k_max - 3.0
        # ymax = k_max
        # xdata, ydata = np.mgrid[xmin:xmax:d, ymin:ymax:d]



        self.F = np.vectorize(lambda x, y, b, f, n: 
                              self.me.cap2_pdf_cdf(b, f, x,
                                                   self.me.k_of_p_tilde(y, k_max),
                                                   self.me.k_of_p_tilde(self._ptilde_B, k_max))[n])

        xdata, ydata = np.mgrid[0.0:(self.me.pars.L0_A*self.me.pars.c0_A):d, 0.01:0.99:d]
        z0data = self.F(xdata, ydata, 0, 0, 0)

        xdata, ydata = np.mgrid[0.0:0.15:d, 0.01:0.99:d]
        z1data = self.F(xdata, ydata, 0, 1, 0)

        xdata, ydata = np.mgrid[0.0:1.0:d, 0.0:1.0:d]

        # x0 = xdata.min()
        # x1 = xdata.max()
        # xdata = (xdata - x0) / (x1 - x0)

        # y0 = ydata.min()
        # y1 = ydata.max()
        # ydata = (ydata - y0) / (y1 - y0)

        z0 = z0data.min()
        z1 = z0data.max()
        z0data = (z0data - z0) / (z1 - z0)

        z0 = z1data.min()
        z1 = z1data.max()
        z1data = (z1data - z0) / (z1 - z0)

        if self.first:
            self.first = False
            self.surfa = self.scenea.mlab.surf(xdata, ydata, z0data, figure=self.scenea.mayavi_scene)#, warp_scale='auto')
            self.surfb = self.sceneb.mlab.surf(xdata, ydata, z1data, figure=self.sceneb.mayavi_scene)#, warp_scale='auto')

            my.axes(self.surfa, nb_labels=5, xlabel="C-ratio c2", ylabel="%-A forecl.", zlabel="Ratio pdf.", figure=self.scenea.mayavi_scene)
            my.axes(self.surfb, nb_labels=5, xlabel="C2", ylabel="%-A forecl.", zlabel="Cap. pdf.", figure=self.sceneb.mayavi_scene)

            # my.clf(figure=self.scenea.mayavi_scene)
            # my.text(0.1, 0.1, "Omega A", figure=self.scenea.mayavi_scene)
            # my.text(0.1, 0.1, "Omega B", figure=self.sceneb.mayavi_scene)
            # my.text(0.1, 0.1, "Omega A+B", figure=self.scenec.mayavi_scene)

            # my.axes(self.surfa, nb_labels=5, xlabel="k_fcl_A", ylabel="k_fcl_B", zlabel="Omega_A", figure=self.scenea.mayavi_scene)
            # my.axes(self.surfb, nb_labels=5, xlabel="k_fcl_A", ylabel="k_fcl_B", zlabel="Omega_B", figure=self.sceneb.mayavi_scene)
            # my.axes(self.surfc, nb_labels=5, xlabel="k_fcl_A", ylabel="k_fcl_B", zlabel="Omega_A+B", figure=self.scenec.mayavi_scene)

            my.sync_camera(self.scenea.mayavi_scene, self.sceneb.mayavi_scene)
            # my.sync_camera(self.scenea.mayavi_scene, self.scenec.mayavi_scene)

        else:
            self.surfa.mlab_source.set(x=xdata, y=ydata, scalars=z0data)
            self.surfb.mlab_source.set(x=xdata, y=ydata, scalars=z1data)

        # Resume rendering
        self.scenea.disable_render = False
        self.sceneb.disable_render = False


        print "."

        # Update done.


    # The layout of the dialog created
    view = View(VSplit(HSplit(Item('scenea', editor=SceneEditor(scene_class=MayaviScene),
                                   height=250, width=300, show_label=False),
                              Item('sceneb', editor=SceneEditor(scene_class=MayaviScene),
                                   height=250, width=300, show_label=False)),
                       #->
                       HSplit(Group('_phi', '_varphi', '_k0',
                                    '_L0_A', '_L0_B', '_rho', '_xi',
                                    '_mu', '_zeta'),
                              Group('_chi',
                                    '_c0', '_c_bar', '_X1', '_ptilde_B', '_mesh_d'),
                              Group('_psi', '_gamma0',
                                    '_gamma1', '_gamma2'))
                       # <-
                       # Group('_', '_phi', '_varphi', '_k0',
                       #       '_L0_A', '_L0_B', '_rho', '_xi',
                       #       '_mu', '_zeta', '_psi', '_chi',
                       #       '_gamma0',
                       #       '_gamma1', '_gamma2', '_q2', '_X1', '_mesh_d'),
                       ),
                resizable=True,
                )

memod = MEModel()
memod.configure_traits()
