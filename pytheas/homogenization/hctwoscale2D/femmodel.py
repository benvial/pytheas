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

    def compute_solution(self, res_list=None):
        res_list = res_list or ["electrostat_scalar_host", "spectral_problem_incl"]
        super().compute_solution(res_list=res_list)

    def compute_modes(self, **kwargs):
        self.analysis = "modal"
        self.compute_solution()
        self.analysis = "direct"

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
        self._print_progress("Retrieving expansion coefficients")
        self.postprocess("postop_coefs_mu")
        filename = self.tmppath("Coefs_mu.txt")
        return femio.load_timetable(filename)

    def postpro_effective_permeability(self):
        """Returns: a function of wavennumber k"""
        # spec = fem.get_spectral_elements()
        kn = self.postpro_eigenvalues()
        norms = self.postpro_norm_eigenvectors()
        coefs = self.postpro_coefs_mu()
        Vincl = self.get_vol_incl()
        coefs_mu = (np.abs(coefs) ** 2) / np.abs((norms)) ** 2

        def mu_eff(k):
            mu_hom = 1
            for i in range(self.neig):
                mu_hom += k ** 2 * coefs_mu[i] / (kn[i] ** 2 - k ** 2)
            return mu_hom

        return mu_eff
