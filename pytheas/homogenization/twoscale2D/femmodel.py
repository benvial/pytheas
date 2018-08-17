import os
import subprocess
import numpy as np
import scipy as sc
from ...tools import femio
from ...basefem import BaseFEM

pi = np.pi

class FemModel(BaseFEM):
    """A class for the two scale convergence homogenization of
    a 2D medium using a finite element model. See the base class :class:`BaseFEM`
    documentation for more info.
    """

    dir_path = os.path.dirname(os.path.abspath(__file__))


    #: flt: caracteristic length of the problem (typically the period)
    l_carac = 1.0


    #: flt: global mesh parameter
    #: `MeshElementSize = l_carac/(parmesh*n)`, `n`: refractive index
    parmesh = 10.

    # opto-geometric parameters  -------------------------------------------
    dx = 1  #: flt: period x
    dy = 1  #: flt: period y
    eps_host = 1-0j
    eps_incl = 1-0j
    dom_des = 1000  #: design domain number (check .geo/.pro files)

    y_flag = False
    save_solution = False

    @property
    def corners_des(self):
        return -self.dx / 2, +self.dx / 2, -self.dy / 2, +self.dy / 2

    def make_param_dict(self):
        param_dict = super().make_param_dict()
        param_dict["y_flag"] = int(self.y_flag)
        param_dict["save_solution"] = int(self.save_solution)
        return param_dict

    def compute_solution(self, **kwargs):
        if self.pattern:
            self.update_epsilon_value()
        self.update_params()
        self.print_progress("Computing solution: homogenization problem")
        argstr = "-petsc_prealloc 1500 -ksp_type preonly \
                 -pc_type lu -pc_factor_mat_solver_package mumps"
        resolution = "electrostat_scalar"
        femio.solve_problem(
            resolution,
            self.path_pro,
            self.path_mesh,
            verbose=self.getdp_verbose,
            path_pos=self.path_pos,
            argstr=argstr,
        )



    def postprocessing(self):
        self.print_progress("Postprocessing")
        subprocess.call(self.ppstr("postop"), shell=True)

    def get_phi(self):
        phi = np.zeros((2, 2), dtype=complex)
        phi[0, 0] = femio.load_table(self.tmp_dir + "/Phixx.txt")
        phi[0, 1] = femio.load_table(self.tmp_dir + "/Phixy.txt")
        phi[1, 0] = femio.load_table(self.tmp_dir + "/Phiyx.txt")
        phi[1, 1] = femio.load_table(self.tmp_dir + "/Phiyy.txt")
        return phi

    def get_vol(self):
        return femio.load_table(self.tmp_dir + "/Vol.txt")

    def get_int_inveps(self):
        Ixx = femio.load_table(self.tmp_dir + "/I_inveps_xx.txt")
        Iyy = femio.load_table(self.tmp_dir + "/I_inveps_yy.txt")
        return Ixx, Iyy

    def postpro_effective_permittivity(self):
        phi = self.get_phi()
        int_inveps_xx, int_inveps_yy = self.get_int_inveps()
        V = self.get_vol()
        epsinv_eff = (np.diag([int_inveps_xx, int_inveps_yy]) + phi) / V
        eps_eff = np.linalg.inv(epsinv_eff)
        return eps_eff

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
        return self.postpro_effective_permittivity()

    def postpro_fields(self, filetype="txt"):
        """ Compute the field maps and output to a file.

            Parameters
            ----------
            filetype : str, default "txt"
                Type of output files. Either "txt" (to be read by the method
                get_field_map in python) or "pos" to be read by gmsh/getdp.

        """
        self.print_progress("Postprocessing fields")
        self.postpro_choice("postop_fields", filetype)
