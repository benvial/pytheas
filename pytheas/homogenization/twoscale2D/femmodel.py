import numpy as np
from ...tools import femio
from ...basefem import *


class TwoScale2D(BaseFEM):
    """A class for the two scale convergence homogenization of
    a 2D medium using a finite element model. See the base class :class:`BaseFEM`
    documentation for more info.
    """

    def __init__(self):
        super().__init__()
        self.dir_path = get_file_path(__file__)

        #: flt: characteristic length of the problem (typically the period)
        self.l_carac = 1.0

        #: flt: global mesh parameter
        #: `MeshElementSize = l_carac/(parmesh*n)`, `n`: refractive index
        self.parmesh = 10.0

        # opto-geometric parameters  -------------------------------------------
        self.dx = 1  #: flt: period x
        self.dy = 1  #: flt: period y
        self.eps_host = 1 - 0j
        self.eps_incl = 1 - 0j
        self.dom_des = 1000  #: design domain number (check .geo/.pro files)

        self.y_flag = False
        self.save_solution = False
        self.aniso = False

    @property
    def corners_des(self):
        return -self.dx / 2, +self.dx / 2, -self.dy / 2, +self.dy / 2

    @property
    def domX_L(self):
        return -self.dx / 2

    @property
    def domX_R(self):
        return self.dx / 2

    @property
    def domY_B(self):
        return -self.dy / 2

    @property
    def domY_T(self):
        # + self.h_pmltop
        return self.dy / 2

    def make_param_dict(self):
        param_dict = super().make_param_dict()
        param_dict["y_flag"] = int(self.y_flag)
        param_dict["save_solution"] = int(self.save_solution)
        param_dict["aniso"] = int(self.aniso)
        return param_dict

    def compute_solution(self, res_list=None):
        res_list = res_list or ["electrostat_scalar", ""]
        super().compute_solution(res_list=res_list)

    def postprocessing(self):
        self.print_progress("Postprocessing")
        self.postprocess("postop")

    def get_phi(self):
        phi = np.zeros((2, 2), dtype=complex)
        phi[0, 0] = femio.load_table(self.tmppath("Phixx.txt"))
        phi[0, 1] = femio.load_table(self.tmppath("Phixy.txt"))
        phi[1, 0] = femio.load_table(self.tmppath("Phiyx.txt"))
        phi[1, 1] = femio.load_table(self.tmppath("Phiyy.txt"))
        return phi

    def get_vol(self):
        return femio.load_table(self.tmppath("Vol.txt"))

    def get_int_inveps(self):
        Ixx = femio.load_table(self.tmppath("I_inveps_xx.txt"))
        Iyy = femio.load_table(self.tmppath("I_inveps_yy.txt"))
        return Ixx, Iyy

    def postpro_effective_permittivity(self):
        phi = self.get_phi()
        int_inveps_xx, int_inveps_yy = self.get_int_inveps()
        if self.python_verbose:
            print("int_inveps_xx = ", int_inveps_xx)
            print("phi = ", phi)
        V = self.get_vol()
        epsinv_eff = (np.diag([int_inveps_xx, int_inveps_yy]) + phi) / V
        # eps_eff = np.linalg.inv(epsinv_eff)
        eps_eff = epsinv_eff.T / np.linalg.det(epsinv_eff)
        return eps_eff

    def get_field_map(self, name):
        field = femio.load_table(self.tmp_dir + "/" + name)
        return np.flipud(field.reshape((self.Niy, self.Nix))).T

    def compute_epsilon_eff(self, postpro_fields=False):
        self.y_flag = False
        self.compute_solution()
        self.postprocessing()
        if postpro_fields:
            self.postpro_fields(filetype="pos")
        self.y_flag = True
        self.compute_solution()
        self.postprocessing()
        if postpro_fields:
            self.postpro_fields(filetype="pos")
        eps_eff = self.postpro_effective_permittivity()
        if self.python_verbose:
            print("#" * 33)
            print("effective permittivity tensor: \n", eps_eff)
        return eps_eff
