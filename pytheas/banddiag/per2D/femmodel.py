# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


import os
import subprocess
import numpy as np
import scipy as sc
from ...tools import femio
from ...basefem import BaseFEM


pi = np.pi


class BandsFEM2D(BaseFEM):
    """A class to compute the band diagramm of a 2D
    bi-periodic medium.
    """

    dir_path = os.path.dirname(os.path.abspath(__file__))

    #: flt: incident wavelength in free space
    lambda0 = 1.0
    lambda_mesh = 1.0
    kx = 0  #: flt: wavevector, x component
    ky = 0  #: flt: wavevector, y component
    eps_des = 1
    analysis = "modal"

    #: str: polarisation (either "TE" or "TM")
    pola = "TE"

    #: flt: shift for eigenvalue search
    lambda0search = 1.0
    #: int: number of eigenvalues
    neig = 6

    #: int: number of x points for postprocessing field maps
    Nix = 51

    #: flt: global mesh parameter
    #: `MeshElementSize = lambda0/(parmesh*n)`, `n`: refractive index
    parmesh = 10.0
    quad_mesh_flag = False
    type_des = "elements"
    y_flag = False
    save_solution = False

    # opto-geometric parameters  -------------------------------------------
    dx = 1  #: flt: period x
    dy = 1  #: flt: period y
    dom_des = 1000  #: design domain number (check .geo/.pro files)

    # postprocessing -------------------------------------------------
    # @property
    # def parmesh_des(self):
    #     return self.parmesh

    @property
    def corners_des(self):
        return -self.dx / 2, +self.dx / 2, -self.dy / 2, +self.dy / 2

    def compute_solution(self, **kwargs):
        if self.pattern:
            self.update_epsilon_value()
        self.update_params()
        self.print_progress("Computing solution")
        argstr = "-slepc -eps_type krylovschur \
                   -st_ksp_type preonly \
                   -st_pc_type lu \
                   -st_pc_factor_mat_solver_package mumps \
                   -eps_max_it 300 \
                   -eps_target 0.00001 \
                   -eps_target_real \
                   -eps_mpd 600 -eps_nev 400"
        resolution = self.pola
        femio.solve_problem(
            resolution,
            self.path_pro,
            self.path_mesh,
            verbose=self.getdp_verbose,
            path_pos=self.path_pos,
            argstr=argstr,
        )

    def postpro_fields(self, filetype="txt"):
        self.print_progress("Postprocessing fields")
        self.postpro_choice("postop_fields", filetype)

    def get_field_map(self, name):
        field = femio.load_table(self.tmp_dir + "/" + name)
        return field.reshape((self.Niy, self.Nix)).T

    def postpro_eigenvalues(self):
        self.print_progress("Retrieving eigenvalues")
        subprocess.call(self.ppstr("postop_eigenvalues_" + self.pola), shell=True)
        filename = self.tmp_dir + "/EV_" + self.pola + ".txt"
        re = np.loadtxt(filename, usecols=[1])
        im = np.loadtxt(filename, usecols=[5])
        return re + 1j * im

    # def postpro_eigenvectors(self, filetype="txt"):
    #     self.print_progress("Retrieving eigenvectors")
    #     self.postpro_choice("postop_eigenvectors_" + self.pola, filetype)
    #     if filetype is "txt":
    #         filename = self.tmp_dir + "/modes_" + self.pola + ".txt"
    #         ev = femio.load_timetable(filename)
    #         return ev.reshape((self.Nix, self.Niy, self.neig))

    def postpro_eigenvectors(self, filetype="txt"):
        self.print_progress("Retrieving eigenvectors")
        self.postpro_choice("postop_eigenvectors_" + self.pola, filetype)
        if filetype is "txt":
            filename = self.tmp_dir + "/modes_" + self.pola + ".txt"
            mode = femio.load_timetable(filename)
            u1 = np.zeros((self.Nix, self.Niy, self.neig), dtype=complex)
            u = mode.reshape((self.Niy, self.Nix, self.neig))
            for imode in range(self.neig):
                u1[:, :, imode] = np.flipud(u[:, :, imode]).T
            return u1
        else:
            return

    def plot_field_and_pattern(
        self,
        fig,
        ax,
        field,
        pattern,
        cmap_div,
        cmap_mat,
        cbar=True,
        vmin=None,
        vmax=None,
    ):

        self.print_progress("Plotting field map")
        x = np.linspace(
            self.nper * self.domX_L, self.nper * self.domX_R, self.nper * self.Nix
        )
        y = np.linspace(self.domY_B, self.domY_T, self.Niy)
        # xx, yy = np.meshgrid(x, y)
        extent = (
            self.nper * self.domX_L,
            self.nper * self.domX_R,
            self.domY_B,
            self.domY_T,
        )

        im1 = ax.imshow(
            field.T,
            interpolation="bilinear",
            cmap=cmap_div,
            vmin=vmin,
            vmax=vmax,
            extent=extent,
        )
        if cbar:
            fig.colorbar(im1, fraction=0.046, pad=0.04)
        ax.imshow(
            pattern.T,
            interpolation="None",
            cmap=cmap_mat,
            alpha=0.33,
            extent=(
                self.nper * self.domX_L,
                self.nper * self.domX_R,
                self.h_layer1,
                self.h_layer1 + self.h_des,
            ),
        )
        ax.imshow(
            field,
            alpha=0.0,
            interpolation="bilinear",
            cmap=cmap_div,
            vmin=vmin,
            vmax=vmax,
            extent=(
                self.nper * self.domX_L,
                self.nper * self.domX_R,
                self.domY_B,
                self.domY_T,
            ),
        )
        ax.set_ylim((self.domY_B, self.domY_T))

    def points_kspace(self, N):

        kx0 = pi / self.dx
        ky0 = pi / self.dy
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

    #
    #
    # def points_kspace(self, N):
    #     Gamma = [0., 0.]
    #     X = [1., 0.]
    #     M = [1., 1.]
    #     ngx = N
    #     Gamma_X = np.array(
    #         [np.linspace(Gamma[0], X[0], ngx), np.linspace(Gamma[1], X[1], ngx)]
    #     )
    #
    #     nxm = N
    #     X_M = np.array([np.linspace(X[0], M[0], nxm), np.linspace(X[1], M[1], nxm)])
    #
    #     X_M = np.delete(X_M, 0, axis=1)
    #
    #     nmg = N
    #     M_Gamma = np.array(
    #         [np.linspace(M[0], Gamma[0], nmg), np.linspace(M[1], Gamma[1], nmg)]
    #     )
    #     M_Gamma = np.delete(M_Gamma, 0, axis=1)
    #     bandsx = np.append(np.append(Gamma_X[0, :], X_M[0, :]), M_Gamma[0, :])
    #     bandsy = np.append(np.append(Gamma_X[1, :], X_M[1, :]), M_Gamma[1, :])
    #
    #     K1 = np.linspace(0, self.dx, ngx)
    #
    #     K2 = np.linspace(self.dx, self.dx + self.dy, nxm)
    #     K2 = np.delete(K2, 0)
    #     K3 = np.linspace(
    #         self.dx + self.dy,
    #         self.dx + self.dy + np.sqrt(self.dx ** 2 + self.dy ** 2),
    #         nmg,
    #     )
    #     K3 = np.delete(K3, 0)
    #     Kplot = np.append(np.append(K1, K2), K3)
    #     K = np.array((bandsx, bandsy)).T
    #     return K, Kplot
