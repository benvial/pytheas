import numpy as np
from ...tools import femio
from ...basefem import *


class TwoScale3D(BaseFEM):
    """A class for the two scale convergence homogenization of
    a 3D medium using a finite element model. See the base class :class:`BaseFEM`
    documentation for more info.
    """

    def __init__(
        self,
        #: flt: characteristic length of the problem (typically the period)
        l_carac=1.0,
        #: flt: global mesh parameter
        #: `MeshElementSize = l_carac/(parmesh*n)`, `n`: refractive index
        parmesh=10.0,
        # opto-geometric parameters  -------------------------------------------
        dx=1,  #: flt: period x
        dy=1,  #: flt: period y
        dz=1,  #: flt: period z
        ax=0.25,  #: flt: ellipsoid principal axis length x
        ay=0.25,  #: flt: ellipsoid principal axis length y
        az=0.25,  #: flt: ellipsoid principal axis length z
        eps_host=1 - 0j,
        eps_incl=1 - 0j,
        dom_des=1002,  #: design domain number (check .geo/.pro files)
        dim=3,  #: dimension of the problem
        save_solution=False,
        type_des="nodes",
        inclusion_flag=False,
        coupling_flag=False,
    ):

        super().__init__()
        self.dir_path = get_file_path(__file__)

        #: flt: characteristic length of the problem (typically the period)
        self.l_carac = l_carac

        #: flt: global mesh parameter
        #: `MeshElementSize = l_carac/(parmesh*n)`, `n`: refractive index
        self.parmesh = parmesh
        # opto-geometric parameters  -------------------------------------------
        self.dx = dx  #: flt: period x
        self.dy = dy  #: flt: period y
        self.dz = dz  #: flt: period z
        self.ax = ax
        self.ay = ay
        self.az = az
        self.eps_host = eps_host
        self.eps_incl = eps_incl
        self.dom_des = dom_des  #: design domain number (check .geo/.pro files)
        self.dim = dim  #: dimension of the problem

        self.save_solution = save_solution

        self.type_des = type_des
        self.lambda_mesh = l_carac
        self.eps_des = eps_host

        self.inclusion_flag = inclusion_flag
        self.coupling_flag = coupling_flag

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

    def _make_param_dict(self):
        param_dict = super()._make_param_dict()
        # param_dict["save_solution"] = int(self.save_solution)
        param_dict["coupling_flag"] = int(self.coupling_flag)
        return param_dict

    def compute_solution(self, **kwargs):
        if self.pattern:
            self.update_epsilon_value()
        self.update_params()
        self._print_progress("Computing solution: homogenization problem")
        argstr = "-petsc_prealloc 1500 -ksp_type preonly \
                 -pc_type lu -pc_factor_mat_solver_type mumps"
        for coord in ["x", "y", "z"]:
            resolution = "annex_" + coord
            self._print_progress("     annex problem " + coord)
            femio.solve_problem(
                resolution,
                self.path_pro,
                self.path_mesh,
                verbose=self.getdp_verbose,
                path_pos=self.path_pos,
                argstr=argstr,
            )
            postop = "postop_" + coord
            self.postprocess(postop)

    def compute_epsilon_eff(self):
        self._print_progress("Postprocessing")
        phi = np.zeros((3, 3), dtype=complex)
        int_eps = []
        for i1, ax1 in enumerate(["x", "y", "z"]):
            int_eps.append(
                femio.load_timetable(self.tmppath("Int_eps_{}{}.txt".format(ax1, ax1)))
            )
            for i2, ax2 in enumerate(["x", "y", "z"]):
                phi[i1, i2] = femio.load_timetable(
                    self.tmppath("Phi" + ax1 + ax2 + ".txt")
                )

        Vcell = self.dx * self.dy * self.dz
        eps_eff = 1 / Vcell * (np.diag(int_eps) - phi)
        # print("int_eps = {}".format(int_eps))
        # print("phi = {}".format(phi))
        if self.python_verbose:
            print("#" * 33)
            print("effective permittivity tensor: \n", eps_eff)
        return eps_eff
