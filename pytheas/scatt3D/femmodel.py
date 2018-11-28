# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


import os
import subprocess
import numpy as np
import scipy as sc
from ..tools import femio
from ..basefem import BaseFEM

from .geom import make_geom

pi = np.pi


class ScattFEM3D(BaseFEM):
    """A class for a finite element model of a 3D scattering object
        using Gmsh_ and GetDP_.

        .. _Gmsh:
            http://gmsh.info/
        .. _GetDP:
            http://getdp.info/
    """

    dir_path = os.path.dirname(os.path.abspath(__file__))

    def __init__(
        self,
        analysis="diffraction",
        A=1,
        lambda0=1,
        theta_deg=0.0,
        phi_deg=0,
        psi_deg=0,
        hx_des=1,  #: flt: thickness x design
        hy_des=1,  #: flt: thickness y design
        hz_des=1,  #: flt: thickness z design
        hx_box=2,  #: flt: thickness x box
        hy_box=2,  #: flt: thickness y box
        hz_box=2,  #: flt: thickness z box
        x_sph=0,  #: flt: sphere center x
        y_sph=0,  #: flt: sphere center y
        z_sph=0.7,  #: flt: sphere center z
        R_sph=0.1,  #: flt: sphere center z
        h_pml=1,  #: flt: thickness PMLs
        a_pml=1,  #: flt: PMLs parameter, real part
        b_pml=1,  #: flt: PMLs parameter, imaginary part
        eps_des=1 - 0 * 1j,  #: flt: permittivity design
        eps_host=1 - 0 * 1j,  #: flt: permittivity host
        eps_sph=1 - 0 * 1j,  #: flt: permittivity sphere
        el_order=1,
    ):
        super().__init__()

        self.analysis = analysis
        self.A = A
        self.lambda0 = lambda0
        self.theta_deg = theta_deg
        self.phi_deg = phi_deg
        self.psi_deg = psi_deg

        # opto-geometric parameters  -------------------------------------------
        #: flt: PMLs parameter, real part
        self.a_pml = a_pml  #: flt: PMLs parameter, real part
        self.b_pml = b_pml  #: flt: PMLs parameter, imaginary part
        self.eps_des = eps_des
        self.eps_host = eps_host
        self.eps_sph = eps_sph
        self.hx_des = hx_des
        self.hy_des = hy_des
        self.hz_des = hz_des
        self.hx_box = hx_box
        self.hy_box = hy_box
        self.hz_box = hz_box
        self.h_pml = h_pml

        self.x_sph = x_sph
        self.y_sph = y_sph
        self.z_sph = z_sph
        self.R_sph = R_sph

        self.el_order = el_order
        self.bg_mesh = False

        # 2  #: design domain number (check .geo/.pro files)
        self.dom_des = 2

        # postprocessing -------------------------------------------------
        #: int: number of x integration points
        #: for postprocessing diffraction efficiencies
        self.ninterv_integ = 60
        #: int: number of z slices points
        #: for postprocessing diffraction efficiencies
        self.nb_slice = 3
        #: flt:  such that `scan_dist  = min(h_sup, hsub)/scan_dist_ratio`
        self.scan_dist_ratio = 5

        self.dim = 3

        self.adjoint = False
        self.recomb_des = True

    @property
    def celltype(self):
        if self.recomb_des:
            make_geom(ext=True)
            return "hexahedron"
        else:
            make_geom(ext=False)
            return "tetra"

    @property
    def theta_0(self):
        return pi / 180.0 * (self.theta_deg)

    @property
    def phi_0(self):
        return pi / 180.0 * (self.phi_deg)

    @property
    def psi_0(self):
        return pi / 180.0 * (self.psi_deg)

    @property
    def corners_des(self):
        return (
            -self.hx_des / 2,
            +self.hx_des / 2,
            -self.hy_des / 2,
            +self.hy_des / 2,
            -self.hz_des / 2,
            +self.hz_des / 2,
        )

    def make_param_dict(self):
        param_dict = super().make_param_dict()
        return param_dict

    def compute_solution(self, **kwargs):
        self.update_params()
        self.print_progress("Computing solution: this may take a while")
        if self.analysis == "diffraction":
            argstr = "-petsc_prealloc 500 -ksp_type preonly \
                     -pc_type lu -pc_factor_mat_solver_package mumps"

            resolution = "helmholtz_vector"
        elif self.analysis == "modal":
            argstr = "-slepc -eps_type krylovschur \
                       -st_ksp_type preonly \
                       -st_pc_type lu \
                       -st_pc_factor_mat_solver_package mumps \
                       -eps_max_it 300 \
                       -eps_target 0.00001 \
                       -eps_target_real \
                       -eps_mpd 600 -eps_nev 400"

            resolution = "helmholtz_vector_modal"
        else:
            raise TypeError(
                "Wrong analysis specified: choose between diffraction and modal"
            )
        argstr += " -cpu"
        femio.solve_problem(
            resolution,
            self.path_pro,
            self.path_mesh,
            verbose=self.getdp_verbose,
            path_pos=self.path_pos,
            argstr=argstr,
        )

    def postpro_absorption(self):
        subprocess.call(self.ppstr("postopQ"), shell=True)
        path = self.tmp_dir + "/Q.txt"
        Q = np.loadtxt(path, skiprows=0, usecols=[1]) + 1j * np.loadtxt(
            path, skiprows=0, usecols=[1]
        )
        return Q.real

    def postpro_epsilon(self):
        self.print_progress("Postprocessing permittivity")
        subprocess.call([self.ppstr("postop_epsilon") + " -order 2"], shell=True)

    #
    def postpro_fields(self, filetype="txt"):
        self.postpro_choice("postop_fields", filetype)

    #
    # def get_field_map(self, name):
    #     field = femio.load_table(self.tmp_dir + "/" + name)
    #     return field.reshape((self.Niy, self.Nix)).T
    #
    def get_objective(self):
        self.print_progress("Retrieving objective")
        if not self.adjoint:
            subprocess.call(self.ppstr("postop_int_objective"), shell=True)
        return femio.load_table(self.tmp_dir + "/objective.txt").real

    def get_adjoint(self):
        self.print_progress("Retrieving adjoint")
        return self.get_qty_vect("adjoint.txt")

    def get_deq_deps(self):
        self.print_progress("Retrieving dEq_deps")
        return self.get_qty_vect("dEq_deps.txt")
