
import numpy as np
from ...tools import femio
from ..twoscale2D.femmodel import TwoScale2D
from ...basefem import get_file_path


class HighContrast2D(TwoScale2D):
    """A class for the two scale convergence homogenization of
    a 2D medium with high contrast using a finite element model. See the base class :class:`BaseFEM`
    documentation for more info.
    """

    def __init__(self):
        super().__init__()
        self.dir_path = get_file_path(__file__)
        self.inclusion_flag = True
        self.inclusion_filename_ = "inclusion.geo"

        self.neig = 10
        self.lambda0search = 100

    def compute_modes(self, **kwargs):
        if self.pattern:
            self.update_epsilon_value()
        self.update_params()
        self.print_progress("Computing solution: spectral problem")
        argstr = "-slepc -eps_type krylovschur \
                   -st_ksp_type preonly \
                   -st_pc_type lu \
                   -st_pc_factor_mat_solver_package mumps \
                   -eps_max_it 300 \
                   -eps_target 0.00001 \
                   -eps_target_real \
                   -eps_mpd 600 -eps_nev 400"
        resolution = "spectral_problem_incl"
        femio.solve_problem(
            resolution,
            self.path_pro,
            self.path_mesh,
            verbose=self.getdp_verbose,
            path_pos=self.path_pos,
            argstr=argstr,
        )

    @property
    def domX_L(self):
        return self.corners_des[0]

    @property
    def domX_R(self):
        return self.corners_des[1]

    @property
    def domY_B(self):
        return self.corners_des[2]

    @property
    def domY_T(self):
        return self.corners_des[3]

    def get_vol_host(self):
        return femio.load_table(self.tmppath("Vol.txt"))

    def get_vol_incl(self):
        self.postprocess("postop_V_incl")
        return femio.load_table(self.tmppath("V_incl.txt"))

    def postpro_effective_permittivity(self):
        phi = self.get_phi()
        int_inveps_xx, int_inveps_yy = self.get_int_inveps()
        V = self.get_vol_host()
        epsinv_eff = (np.diag([int_inveps_xx, int_inveps_yy]) + phi) / V
        eps_eff = np.linalg.inv(epsinv_eff)
        return eps_eff

    def postpro_coefs_mu(self):
        self.print_progress("Retrieving expansion coefficients")
        self.postprocess("postop_coefs_mu")
        filename = self.tmppath("Coefs_mu.txt")
        return femio.load_timetable(filename)
