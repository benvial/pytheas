# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


import numpy as np

# from ..tools import femio
from ..basefem import *
from .geom import make_geom


class Scatt3D(BaseFEM):
    """A class for a finite element model of a 3D scattering object
        using Gmsh_ and GetDP_.

        .. _Gmsh:
            http://gmsh.info/
        .. _GetDP:
            http://getdp.info/
    """

    def __init__(
        self,
        analysis="direct",
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
        self.dir_path = get_file_path(__file__)

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
        return np.pi / 180.0 * (self.theta_deg)

    @property
    def phi_0(self):
        return np.pi / 180.0 * (self.phi_deg)

    @property
    def psi_0(self):
        return np.pi / 180.0 * (self.psi_deg)

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

    def _make_param_dict(self):
        param_dict = super()._make_param_dict()
        return param_dict

    def compute_solution(self, **kwargs):
        res_list = ["helmholtz_vector", "helmholtz_vector_modal"]
        return super().compute_solution(res_list=res_list)

    def postpro_epsilon(self):
        self._print_progress("Postprocessing permittivity")
        self.postprocess("postop_epsilon" + " -order 2")
