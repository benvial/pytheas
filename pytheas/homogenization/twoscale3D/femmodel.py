import os
import subprocess
import numpy as np
import scipy as sc
from ...tools import femio
from ...basefem import BaseFEM

pi = np.pi


class FemModel(BaseFEM):
    """A class for the two scale convergence homogenization of
    a 3D medium using a finite element model. See the base class :class:`BaseFEM`
    documentation for more info.
    """

    dir_path = os.path.dirname(os.path.abspath(__file__))

    def __init__(
        self,
        #: flt: caracteristic length of the problem (typically the period)
        l_carac=1.0,
        #: flt: global mesh parameter
        #: `MeshElementSize = l_carac/(parmesh*n)`, `n`: refractive index
        parmesh=10.0,
        # opto-geometric parameters  -------------------------------------------
        dx=1,  #: flt: period x
        dy=1,  #: flt: period y
        dz=1,  #: flt: period z
        eps_host=1 - 0j,
        eps_incl=1 - 0j,
        dom_des=1002,  #: design domain number (check .geo/.pro files)
        dim=3,  #: dimension of the problem
        y_flag=False,
        save_solution=False,
        type_des="nodes",
    ):

        super().__init__()

        #: flt: caracteristic length of the problem (typically the period)
        self.l_carac = l_carac

        #: flt: global mesh parameter
        #: `MeshElementSize = l_carac/(parmesh*n)`, `n`: refractive index
        self.parmesh = parmesh
        # opto-geometric parameters  -------------------------------------------
        self.dx = dx  #: flt: period x
        self.dy = dx  #: flt: period y
        self.dy = dx  #: flt: period z
        self.eps_host = eps_host
        self.eps_incl = eps_incl
        self.dom_des = dom_des  #: design domain number (check .geo/.pro files)
        self.dim = dim  #: dimension of the problem

        self.y_flag = y_flag
        self.save_solution = save_solution

        self.type_des = type_des
        self.lambda_mesh = l_carac
        self.eps_des = eps_host

    celltype = "tetra"

    @property
    def corners_des(self):
        return (
            -self.dx / 2,
            +self.dx / 2,
            -self.dy / 2,
            +self.dy / 2,
            -self.dz / 2,
            +self.dz / 2,
        )

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
        for coord in ["x", "y", "z"]:
            resolution = "electrostat_scalar_" + coord
            femio.solve_problem(
                resolution,
                self.path_pro,
                self.path_mesh,
                verbose=self.getdp_verbose,
                path_pos=self.path_pos,
                argstr=argstr,
            )
            postop = "postop_" + coord
            subprocess.call(self.ppstr(postop), shell=True)

    def compute_epsilon_eff(self):
        self.print_progress("Postprocessing")
        phi = np.zeros((3, 3), dtype=complex)
        i1 = 0
        for ax1 in ("x", "y", "z"):
            i2 = 0
            for ax2 in ("x", "y", "z"):
                phi[i1, i2] = femio.load_timetable(
                    self.tmp_dir + "/Phi" + ax1 + ax2 + ".txt"
                )
                i2 += 1
            i1 += 1
        int_eps = femio.load_timetable(self.tmp_dir + "/Inteps.txt")
        Vcell = self.dx * self.dy * self.dz
        eps_eff = 1 / Vcell * (int_eps * np.eye(3) - phi)
        if self.python_verbose:
            print("#" * 33)
            print("effective permittivity tensor: \n", eps_eff)
        return eps_eff

    #
    # def postpro_fields(self, filetype="txt"):
    #     """ Compute the field maps and output to a file.
    #
    #         Parameters
    #         ----------
    #         filetype : str, default "txt"
    #             Type of output files. Either "txt" (to be read by the method
    #             get_field_map in python) or "pos" to be read by gmsh/getdp.
    #
    #     """
    #     self.print_progress("Postprocessing fields")
    #     self.postpro_choice("postop_fields", filetype)
