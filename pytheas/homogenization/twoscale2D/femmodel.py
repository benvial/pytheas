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
        self.interp = True
        self.N_pts_circ = 100
        self.R_circ = 1

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

    def _make_param_dict(self):
        param_dict = super()._make_param_dict()
        param_dict["y_flag"] = int(self.y_flag)
        param_dict["save_solution"] = int(self.save_solution)
        param_dict["aniso"] = int(self.aniso)
        param_dict["interp"] = int(self.interp)
        return param_dict

    def compute_solution(self, res_list=None):
        res_list = res_list or ["electrostat_scalar", ""]
        super().compute_solution(res_list=res_list)

    def postprocessing(self):
        self._print_progress("Postprocessing")
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

    def get_solution(self):
        return femio.load_table(self.tmppath("u.txt"))

    def get_int_inveps(self):
        I = np.zeros((2, 2), dtype=complex)
        I[0, 0] = femio.load_table(self.tmppath("I_inveps_xx.txt"))
        I[1, 1] = femio.load_table(self.tmppath("I_inveps_yy.txt"))
        I[1, 0] = femio.load_table(self.tmppath("I_inveps_yx.txt"))
        I[0, 1] = femio.load_table(self.tmppath("I_inveps_xy.txt"))
        return I

    def postpro_effective_permittivity(self):
        V = self.get_vol()
        int_inveps = self.get_int_inveps() / V
        phi = self.get_phi() / V
        if self.python_verbose:
            print("int_inveps = ", int_inveps)
            print("phi = ", phi)
        epsinv_eff = int_inveps + phi
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

    def get_deq_deps(self):
        return self._get_qty("dEq_deps_x.txt"), self._get_qty("dEq_deps_y.txt")

    def postpro_circ(self):
        for s in ["sol", "vx", "vy"]:
            t = self.tmppath("{}_circ.txt".format(s))
            try:
                os.remove(t)
            except:
                pass
        self.postprocess("postop_fields_circle")

    def get_sol_circ(self):
        t = self.tmppath("sol_circ.txt")
        return femio.load_table(t)

    def get_gradsol_circ(self):
        vx = femio.load_table(self.tmppath("vx_circ.txt"))
        vy = femio.load_table(self.tmppath("vy_circ.txt"))
        return vx, vy

    #
    # def get_laplacian_psi(self, interp_method="nearest"):
    #     deq_deps = self._get_qty("dx_psi.txt"), self._get_qty("dy_psi.txt")
    #     x_grid, y_grid = self.grid
    #
    #     deq_deps_x, deq_deps_y = deq_deps
    #     deq_deps_x = self.mesh2grid(deq_deps_x, interp_method=interp_method)
    #     deq_deps_x_x = np.gradient(deq_deps_x.T)[0] / np.gradient(x_grid)[0]
    #
    #     deq_deps_y = self.mesh2grid(deq_deps_y, interp_method=interp_method)
    #     deq_deps_y_y = np.gradient(deq_deps_y.T)[1] / np.gradient(y_grid)[0]
    #
    #     deq_deps = deq_deps_x_x.T + deq_deps_y_y.T
    #
    #     deq_deps_re = self.grid2mesh(deq_deps.real)
    #     deq_deps_im = self.grid2mesh(deq_deps.imag)
    #     deq_deps = deq_deps_re + 1j * deq_deps_im
    #     return deq_deps
