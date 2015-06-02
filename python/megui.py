

import numpy as np
import mayavi.mlab as my

from traits.api import HasTraits, Range, Instance, \
        on_trait_change
from traitsui.api import View, Item, Group, HSplit, VSplit

from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import MayaviScene, SceneEditor, \
                MlabSceneModel


# Load the library containing an interface to our C solver
import pybnet


# **********************************************************************
# INITIALIZATION


class MEModel(HasTraits):
    _rho = Range(0.05, 0.95, 0.3, )
    _gamma = Range(0.0, 1.0, 0.4, )
    _psi = Range(0.0, 10.0, 1.0, )
    _chi = Range(0.0, 1.0, 0.0, )
    _theta = Range(0.0, 1.0, 0.2, )
    _ca = Range(0.01, 0.2, 0.06, )
    _cb = Range(0.01, 0.2, 0.06, )
    _La = Range(0.1, 10.0, 4.5, )
    _Lb = Range(0.1, 10.0, 3.0, )
    _X = Range(-3.0, 1.0, -1.5, )

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
        self.bn = pybnet.Bnet(LGD = 0.2,
                              PD = 0.05,
                              rho = self._rho,
                              gamma = self._gamma,
                              psi = self._psi,
                              q0 = 0.1,
                              q1 = 0.9,
                              chi = self._chi,
                              theta = self._theta)

        # Calibrate phi, and phi2 relative to phi

    @on_trait_change('scenec.activated')
    def populate_scenes(self):
        self.update_plot()

    @on_trait_change('_rho, _gamma, _psi, _chi, _theta, _ca, _cb, _La, _Lb, _X')
    def update_plot(self):
        
        self.bn.pars.rho = self._rho
        self.bn.pars.gamma = self._gamma
        self.bn.pars.psi = self._psi
        self.bn.pars.chi = self._chi
        self.bn.pars.theta = self._theta

        phi_calib = self.bn.calibrate_phi()
        print("phi_calib = {}".format(phi_calib))

        self.bn.pars.phi = phi_calib
        self.bn.pars.phi2 = 2*self.bn.pars.phi

        # print("Update: rho -> {}, gamma -> {}, psi -> {}, chi -> {}, theta -> {}".format(self._rho, self._gamma, self._psi, self._chi, self._theta))

        kzero = self.bn.comp_kzero()

        ka0 = kzero - 2.5
        ka1 = kzero
        kb0 = kzero - 2.5
        kb1 = kzero

        kadata, kbdata = np.mgrid[ka0:ka1:10j, kb0:kb1:10j]

        self.K = np.vectorize(lambda ka, kb: 
                              self.bn.comp_kimp((self.bn.comp_epsstar(ka, self._X),
                                                 self.bn.comp_epsstar(kb, self._X)),
                                                (self._La, self._Lb), self._X))
        self.WA = np.vectorize(lambda ka, kb, n: 
                               self.bn.comp_Omegaplus(self._ca, self._X, 
                                                      ka, self.K(ka, kb))[n])
        self.WB = np.vectorize(lambda ka, kb, n: 
                               self.bn.comp_Omegaplus(self._cb, self._X, 
                                                      kb, self.K(ka, kb))[n])

        # print("Update ca -> {}, cb -> {}".format(self._ca, self._cb))

        WAdata = self.WA(kadata, kbdata, 0)
        WBdata = self.WB(kadata, kbdata, 0)
        WSdata = WAdata + WBdata

        z0 = WAdata.min()
        z1 = WAdata.max()
        WAdata = (WAdata - z0) / (z1 - z0)

        z0 = WBdata.min()
        z1 = WBdata.max()
        WBdata = (WBdata - z0) / (z1 - z0)

        z0 = WSdata.min()
        z1 = WSdata.max()
        WSdata = (WSdata - z0) / (z1 - z0)


        if self.first:
            self.first = False
            self.surfa = self.scenea.mlab.surf(kadata, kbdata, WAdata, figure=self.scenea.mayavi_scene)#, warp_scale='auto')
            self.surfb = self.sceneb.mlab.surf(kadata, kbdata, WBdata, figure=self.sceneb.mayavi_scene)#, warp_scale='auto')
            self.surfc = self.sceneb.mlab.surf(kadata, kbdata, WSdata, figure=self.scenec.mayavi_scene)#, warp_scale='auto')
            my.axes(self.surfa, nb_labels=5, xlabel="k_A", ylabel="k_B", zlabel="W_A", figure=self.scenea.mayavi_scene)
            my.axes(self.surfb, nb_labels=5, xlabel="k_A", ylabel="k_B", zlabel="W_B", figure=self.sceneb.mayavi_scene)
            my.axes(self.surfc, nb_labels=5, xlabel="k_A", ylabel="k_B", zlabel="W_(A+B)", figure=self.scenec.mayavi_scene)
            # my.title("Scaled Omega A", figure=self.scenea.mayavi_scene)
            # my.title("Scaled Omega B", figure=self.sceneb.mayavi_scene)
            # my.title("Scaled Omega A + Omega B", figure=self.scenec.mayavi_scene)
            my.sync_camera(self.scenea.mayavi_scene, self.sceneb.mayavi_scene)
            my.sync_camera(self.scenea.mayavi_scene, self.scenec.mayavi_scene)

        else:
            self.scenea.disable_render = True
            self.sceneb.disable_render = True
            self.scenec.disable_render = True
            self.surfa.mlab_source.set(x=kadata, y=kbdata, scalars=WAdata)
            self.surfb.mlab_source.set(x=kadata, y=kbdata, scalars=WBdata)  
            self.surfc.mlab_source.set(x=kadata, y=kbdata, scalars=WSdata) 
            self.scenea.disable_render = False
            self.sceneb.disable_render = False
            self.scenec.disable_render = False


    # The layout of the dialog created
    view = View(VSplit(HSplit(Item('scenea', editor=SceneEditor(scene_class=MayaviScene),
                                   height=250, width=300, show_label=False),
                              Item('sceneb', editor=SceneEditor(scene_class=MayaviScene),
                                   height=250, width=300, show_label=False),
                              Item('scenec', editor=SceneEditor(scene_class=MayaviScene),
                                   height=250, width=300, show_label=False),
                              ),
                       Group('_', '_rho', '_gamma', '_psi', '_chi', '_theta', '_X', '_ca', '_cb', '_La', '_Lb',
                             ),
                       ),
                resizable=True,
                )



memod = MEModel()
memod.configure_traits()
