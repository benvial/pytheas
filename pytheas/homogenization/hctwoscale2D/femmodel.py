import os
import subprocess
import numpy as np
import scipy as sc
from ...tools import femio
from ..twoscale2D.femmodel import FemModel as FemModel_

pi = np.pi


class FemModel(FemModel_):
    """A class for the two scale convergence homogenization of
    a 2D medium with high contrast using a finite element model. See the base class :class:`BaseFEM`
    documentation for more info.
    """

    dir_path = os.path.dirname(os.path.abspath(__file__))

    def __init__(self):
        super().__init__()
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
        return femio.load_table(self.tmp_dir + "/Vol.txt")

    def get_vol_incl(self):
        subprocess.call(self.ppcmd("postop_V_incl"))
        return femio.load_table(self.tmp_dir + "/V_incl.txt")

    def postpro_effective_permittivity(self):
        phi = self.get_phi()
        int_inveps_xx, int_inveps_yy = self.get_int_inveps()
        V = self.get_vol_host()
        epsinv_eff = (np.diag([int_inveps_xx, int_inveps_yy]) + phi) / V
        eps_eff = np.linalg.inv(epsinv_eff)
        return eps_eff

    def postpro_eigenvalues(self):
        self.print_progress("Retrieving eigenvalues")
        subprocess.call(self.ppcmd("postop_eigenvalues"))
        filename = self.tmp_dir + "/EigenValues.txt"
        re = np.loadtxt(filename, usecols=[1])
        im = np.loadtxt(filename, usecols=[5])
        return re + 1j * im

    def get_spectral_elements(self):
        eigval = self.postpro_eigenvalues()
        eigvect = self.postpro_eigenvectors()
        isort = np.argsort(eigval)
        eigval = eigval[isort]
        eigvect = eigvect[:, :, (isort)]
        nms = self.postpro_norm_eigenvectors()
        return eigval, eigvect / nms[isort]

    def postpro_eigenvectors(self, filetype="txt"):
        self.print_progress("Retrieving eigenvectors")
        self.postpro_choice("postop_eigenvectors", filetype)
        if filetype is "txt":
            mode = femio.load_timetable(self.tmp_dir + "/EigenVectors.txt")
            u1 = np.zeros((self.Nix, self.Niy, self.neig), dtype=complex)
            u = mode.reshape((self.Niy, self.Nix, self.neig))
            for imode in range(self.neig):
                u1[:, :, imode] = np.flipud(u[:, :, imode]).T
            return u1
        else:
            return

    def postpro_norm_eigenvectors(self):
        self.print_progress("Retrieving eigenvector norms")
        subprocess.call(self.ppcmd("postop_norm_eigenvectors"))
        filename = self.tmp_dir + "/NormsEigenVectors.txt"
        return np.sqrt(femio.load_timetable(filename))

    def postpro_coefs_mu(self):
        self.print_progress("Retrieving expansion coefficients")
        subprocess.call(self.ppcmd("postop_coefs_mu"))
        filename = self.tmp_dir + "/Coefs_mu.txt"
        return femio.load_timetable(filename)
