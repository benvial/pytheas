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

    def __init__(self):
        super().__init__()

        #: flt: caracteristic length of the problem (typically the period)
        self.l_carac = 1.0

        self.inclusion_flag=True

        #: flt: global mesh parameter
        #: `MeshElementSize = l_carac/(parmesh*n)`, `n`: refractive index
        self.parmesh = 10.

        # opto-geometric parameters  -------------------------------------------
        self.dx = 1  #: flt: period x
        self.dy = 1  #: flt: period y
        self.eps_host = 1 - 0j
        self.eps_incl = 1 - 0j
        self.dom_des = 1000  #: design domain number (check .geo/.pro files)

        self.y_flag = False
        self.save_solution = False
        self.inclusion_filename_ = "inclusion.geo"


        self.neig = 10
        self.lambda0search = 100

    @property
    def corners_des(self):
        return -self.dx / 2, +self.dx / 2, -self.dy / 2, +self.dy / 2

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
        resolution = "electrostat_scalar_host"
        femio.solve_problem(
            resolution,
            self.path_pro,
            self.path_mesh,
            verbose=self.getdp_verbose,
            path_pos=self.path_pos,
            argstr=argstr,
        )

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

    def get_vol_host(self):
        return femio.load_table(self.tmp_dir + "/Vol.txt")

    def get_vol_incl(self):
        subprocess.call(self.ppstr("postop_V_incl"), shell=True)
        return femio.load_table(self.tmp_dir + "/V_incl.txt")

    def get_int_inveps(self):
        Ixx = femio.load_table(self.tmp_dir + "/I_inveps_xx.txt")
        Iyy = femio.load_table(self.tmp_dir + "/I_inveps_yy.txt")
        return Ixx, Iyy

    def postpro_effective_permittivity(self):
        phi = self.get_phi()
        int_inveps_xx, int_inveps_yy = self.get_int_inveps()
        V = self.get_vol_host()
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


    def postpro_eigenvalues(self):
        self.print_progress("Retrieving eigenvalues")
        subprocess.call(self.ppstr("postop_eigenvalues"), shell=True)
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
        return eigval, eigvect/nms[isort]

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
        subprocess.call(self.ppstr("postop_norm_eigenvectors"), shell=True)
        filename = self.tmp_dir + "/NormsEigenVectors.txt"
        return np.sqrt(femio.load_timetable(filename))



    def postpro_coefs_mu(self):
        self.print_progress("Retrieving expansion coefficients")
        subprocess.call(self.ppstr("postop_coefs_mu"), shell=True)
        filename = self.tmp_dir + "/Coefs_mu.txt"
        re = np.loadtxt(filename, usecols=[1])
        im = np.loadtxt(filename, usecols=[5])
        return femio.load_timetable(filename)
