# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT

"""
Finite Element model computing the band diagramm of a 2D
bi-periodic medium

"""

import numpy as np
from ...tools import femio
from ...basefem import *


class BandDiag2D(BaseFEM):
    """A class to compute the band diagramm of a 2D
    bi-periodic medium.
    """

    def __init__(self, ID="sim"):
        super().__init__()
        self.dir_path = get_file_path(__file__)

        #: flt: incident wavelength in free space
        self.lambda0 = 1.0
        self.lambda_mesh = 1.0
        self.kx = 0  #: flt: wavevector, x component
        self.ky = 0  #: flt: wavevector, y component
        self.eps_des = 1
        self.analysis = "modal"

        #: str: polarisation (either "TE" or "TM")
        self.pola = "TE"

        #: flt: shift for eigenvalue search
        self.lambda0search = 1.0
        #: int: number of eigenvalues
        self.neig = 6

        #: int: number of x points for postprocessing field maps
        self.Nix = 51

        #: flt: global mesh parameter
        #: `MeshElementSize = lambda0/(parmesh*n)`, `n`: refractive index
        self.parmesh = 10.0
        self.quad_mesh_flag = False
        self.type_des = "elements"
        self.y_flag = False
        self.save_solution = False

        # opto-geometric parameters  -------------------------------------------
        self.dx = 1  #: flt: period x
        self.dy = 1  #: flt: period y
        self.dom_des = 1000  #: design domain number (check .geo/.pro files)

    # postprocessing -------------------------------------------------
    # @property
    # def parmesh_des(self):
    #     return self.parmesh

    @property
    def corners_des(self):
        return -self.dx / 2, +self.dx / 2, -self.dy / 2, +self.dy / 2

    def compute_solution(self):
        super().compute_solution(res_list=[None, self.pola])

    def get_field_map(self, name):
        field = femio.load_table(self.tmppath(name))
        return field.reshape((self.Niy, self.Nix)).T

    def postpro_eigenvalues(self):
        eig_file = "EV_" + self.pola + ".txt"
        postop = "postop_eigenvalues_" + self.pola
        return super().postpro_eigenvalues(postop=postop, eig_file=eig_file)

    def postpro_eigenvectors(self, filetype="txt"):
        postop = "postop_eigenvectors_" + self.pola
        eig_file = "modes_" + self.pola + ".txt"
        return super().postpro_eigenvectors(
            filetype=filetype, postop=postop, eig_file=eig_file
        )

    def points_kspace(self, N):

        kx0 = np.pi / self.dx
        ky0 = np.pi / self.dy
        Gamma = [0.0, 0.0]
        X = [kx0, 0.0]
        M = [kx0, ky0]
        ngx = N
        Gamma_X = np.array(
            [np.linspace(Gamma[0], X[0], ngx), np.linspace(Gamma[1], X[1], ngx)]
        )

        nxm = N
        X_M = np.array([np.linspace(X[0], M[0], nxm), np.linspace(X[1], M[1], nxm)])

        X_M = np.delete(X_M, 0, axis=1)

        nmg = N
        M_Gamma = np.array(
            [np.linspace(M[0], Gamma[0], nmg), np.linspace(M[1], Gamma[1], nmg)]
        )
        M_Gamma = np.delete(M_Gamma, 0, axis=1)
        bandsx = np.append(np.append(Gamma_X[0, :], X_M[0, :]), M_Gamma[0, :])
        bandsy = np.append(np.append(Gamma_X[1, :], X_M[1, :]), M_Gamma[1, :])

        K1 = np.linspace(0, kx0, ngx)
        K2 = kx0 + np.linspace(0, ky0, nxm)
        K2 = np.delete(K2, 0)
        K3 = kx0 + ky0 + np.linspace(0, np.sqrt(kx0 ** 2 + ky0 ** 2), nmg)
        K3 = np.delete(K3, 0)
        Kplot = np.append(np.append(K1, K2), K3)
        K = np.array((bandsx, bandsy)).T
        return K, Kplot
