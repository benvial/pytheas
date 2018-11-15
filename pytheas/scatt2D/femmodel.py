# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT

"""
Finite Element model of 2D media

"""

import os
import subprocess
import numpy as np
import scipy as sc
from ..tools import femio
from ..basefem import BaseFEM


pi = np.pi


class FemModel(BaseFEM):
    """A class for a finite element model of a 2D medium"""

    dir_path = os.path.dirname(os.path.abspath(__file__))

    def __init__(self):
        super().__init__()

        #: str: analysys type (either "diffraction" or "modal")
        self.analysis = "diffraction"
        #: str: polarisation of the incident plane wave
        #: (either "TE" or "TM")
        self.pola = "TE"

        self.nb_incl = 1

        # Incident plane wave parameters
        #: flt: incident plane wave amplitude
        self.A = 1.0
        #: flt: incident plane wave wavelength in free space
        self.lambda0 = 1.0
        #: flt: wavelength to use for meshing
        self.lambda_mesh = 1.0
        #: flt : incident plane wave angle (in degrees).
        #: Light comes from the top
        #: (travels along -y if normal incidence, `theta_deg=0` is set)
        self.theta_deg = 0.

        #: line source position
        self.ls_flag = False
        self.xs = 0
        self.ys = 0

        #: beam?
        self.beam_flag = False
        self.waist = 1.5

        # opto-geometric parameters  -------------------------------------------
        self.h_pml = 2.  #: flt: thickness pml
        self.hx_des = 2.  #: flt: x - thickness scattering box (design)
        self.hy_des = 2.  #: flt: y - thickness scattering box
        self.a_pml = 1  #: flt: PMLs parameter, real part
        self.b_pml = 1  #: flt: PMLs parameter, imaginary part
        self.eps_host = 1  #: flt: permittivity host
        self.eps_sub = 1  #: flt: permittivity substrate
        self.eps_des = 1  #: flt: permittivity scattering box
        self.eps_incl = 1  #: flt: permittivity inclusion
        self.dom_des = 5  #: design domain number (check .geo/.pro files)

        # postprocessing -------------------------------------------------
        #: coords of point for PostProcessing
        self.xpp, self.ypp = 0, 0
        #: int: number of theta points for computing the angular dependance
        #: of the modal coupling coefficients
        self.Ni_theta = 45
        self.M_fs = 3
        self.Nibox_x = 500  #: int: number of x interpolation points on the design box
        self.Nibox_y = 500  #: int: number of y interpolation points on the design box
        self.Nin2f_x = (
            500
        )  #: int: number of x interpolation points for near to far field calculations
        self.Nin2f_y = (
            500
        )  #: int: number of y interpolation points for near to far field calculations
        self.rat_n2f = 0.95

        #: int: number of y slices points
        #: for postprocessing diffraction efficiencies
        self.nb_slice = 20
        #: flt:  such that `scan_dist  = min(h_sup, hsub)/scan_dist_ratio`
        self.scan_dist_ratio = 5

        self.space2pml_L, self.space2pml_R = 1.0, 1.0
        self.space2pml_T, self.space2pml_B = 1.0, 1.0

        #: int: number of x points for postprocessing field maps
        self.Nix = 100
        self.Niy = 100
        self.nper = 3
        # modal analysis -------------------------------------------------
        #: int: number of eigenvalues searched for in modal analysis
        self.neig = 6
        #: flt: wavelength around which to search eigenvalues
        self.lambda0search = 1

        # topology optimization   ------------------------
        self.adjoint = False
        self.target_flag = False
        self.r_target = 0.01
        self.x_target = 0
        self.y_target = 0

    @property
    def corners_des(self):
        return -self.hx_des / 2, +self.hx_des / 2, -self.hy_des / 2, +self.hy_des / 2

    @property
    def domX_L(self):
        return -self.hx_des / 2 - self.space2pml_L

    @property
    def domX_R(self):
        return self.hx_des / 2 + self.space2pml_R

    @property
    def domY_B(self):
        return -self.hy_des / 2 - self.space2pml_B

    @property
    def domY_T(self):
        return self.hy_des / 2 + self.space2pml_T

    @property
    def Xn2f_L(self):
        return self.domX_L * self.rat_n2f

    @property
    def Xn2f_R(self):
        return self.domX_R * self.rat_n2f

    @property
    def Yn2f_B(self):
        return self.domY_B * self.rat_n2f

    @property
    def Yn2f_T(self):
        return self.domY_T * self.rat_n2f

    @property
    def theta(self):
        return pi / 180. * (self.theta_deg)

    @property
    def omega0(self):
        return 2. * pi * self.cel / self.lambda0

    def make_param_dict(self):
        param_dict = super().make_param_dict()
        param_dict["ls_flag"] = int(self.ls_flag)
        param_dict["beam_flag"] = int(self.beam_flag)
        param_dict["target_flag"] = int(self.target_flag)
        return param_dict

    def postpro_fields_box(self):
        path_T = self.tmp_dir + "/field_box_T.out"
        path_L = self.tmp_dir + "/field_box_L.out"
        path_R = self.tmp_dir + "/field_box_R.out"
        path_B = self.tmp_dir + "/field_box_B.out"
        subprocess.call(self.ppstr("postop_fields_box"), shell=True)
        u_T = femio.load_timetable(path_T)
        u_R = femio.load_timetable(path_R)
        u_B = femio.load_timetable(path_B)
        u_L = femio.load_timetable(path_L)
        if self.analysis == "modal":
            u_T = u_T.reshape((self.Nibox_x, self.neig))
            u_B = u_B.reshape((self.Nibox_x, self.neig))
            u_L = u_L.reshape((self.Nibox_y, self.neig))
            u_R = u_R.reshape((self.Nibox_y, self.neig))
        return u_T, u_B, u_L, u_R

    def postpro_fields_n2f(self):
        subprocess.call(self.ppstr("postop_fields_n2f"), shell=True)
        u_out, vx_out, vy_out = {}, {}, {}
        for i in ["T", "L", "R", "B"]:
            u = femio.load_timetable(self.tmp_dir + "/field_n2f_{0}.out".format(i))
            vx = femio.load_timetable(
                self.tmp_dir + "/field_dual_x_n2f_{0}.out".format(i)
            )
            vy = femio.load_timetable(
                self.tmp_dir + "/field_dual_y_n2f_{0}.out".format(i)
            )
            if self.analysis == "modal":
                if (i == "T") or (i == "B"):
                    n = self.Nin2f_x
                else:
                    n = self.Nin2f_y
                u = u.reshape((n, self.neig))
                vx = vx.reshape((n, self.neig))
                vy = vy.reshape((n, self.neig))
            u_out[i] = u
            vx_out[i] = vx
            vy_out[i] = vy
        return u_out, vx_out, vy_out

    def normalized_scs(self, ff, theta):

        xi = np.linspace(self.Xn2f_L, self.Xn2f_R, self.Nin2f_x)
        yi = np.linspace(self.Yn2f_B, self.Yn2f_T, self.Nin2f_y)
        nu0 = np.sqrt(self.mu0 / self.epsilon0)
        k0 = 2 * pi / self.lambda0


        x = {"T": xi, "L": self.Xn2f_L, "R": self.Xn2f_R, "B": xi}
        y = {"T": self.Yn2f_T, "L": yi, "R": yi, "B": self.Yn2f_B}
        nx = {"T": 0, "L": -1, "R": 1, "B": 0}
        ny = {"T": 1, "L": 0, "R": 0, "B": -1}
        Ez, Hx, Hy = ff

        Itheta = []
        for t in theta:
            s = np.sin(t)
            c = -np.cos(t)

            I = 0
            for loc in ["T", "L", "R", "B"]:
                if (loc == "T") or (loc == "B"):
                    l = xi
                else:
                    l = yi
                expo = np.exp(-1j * k0 * (x[loc] * s + y[loc] * c))
                J = (
                    nu0 * (nx[loc] * Hy[loc] - ny[loc] * Hx[loc])
                    - (ny[loc] * c + nx[loc] * s) * Ez[loc]
                ) * expo
                I += np.trapz(J, l)
            Itheta.append(I)
        return np.abs(Itheta) ** 2 / (8 * pi)

    def postpro_fields(self, filetype="txt"):
        self.print_progress("Postprocessing fields")
        self.postpro_choice("postop_fields", filetype)

    def get_field_map(self, name):
        field = femio.load_table(self.tmp_dir + "/" + name)
        return np.flipud(field.reshape((self.Niy, self.Nix)).T)

    def get_field_point(self):
        subprocess.call(self.ppstr("postop_field_on_point"), shell=True)
        u_tot = femio.load_table(self.tmp_dir + "/" + "u_tot_point.txt")
        u_i = femio.load_table(self.tmp_dir + "/" + "u_i_point.txt")
        u = femio.load_table(self.tmp_dir + "/" + "u_point.txt")
        return u, u_tot, u_i

    def get_objective(self):
        self.print_progress("Retrieving objective")
        if not self.adjoint:
            subprocess.call(self.ppstr("postop_int_objective"), shell=True)
        return femio.load_table(self.tmp_dir + "/objective.txt").real

    def postpro_eigenvalues(self):
        self.print_progress("Retrieving eigenvalues")
        subprocess.call(self.ppstr("postop_eigenvalues"), shell=True)
        filename = self.tmp_dir + "/EigenValues.txt"
        re = np.loadtxt(filename, usecols=[1])
        im = np.loadtxt(filename, usecols=[5])
        return re + 1j * im

    def postpro_norm_eigenvectors(self):
        self.print_progress("Retrieving eigenvector norms")
        subprocess.call(self.ppstr("postop_norm_eigenvectors"), shell=True)
        filename = self.tmp_dir + "/NormsEigenVectors.txt"
        return femio.load_timetable(filename)

    def postpro_coupling_angle(self):
        self.print_progress("Angular sweep for coupling coeffs")
        subprocess.call(self.ppstr("postop_coupling_coeffs_angle"), shell=True)
        filename = self.tmp_dir + "/coupling_coeffs.txt"
        tmp = femio.load_timetable(filename)
        return tmp.reshape((self.Ni_theta, self.neig))

    def postpro_fourrier_coefs_angle(self):
        self.print_progress("Fourrier coefficients for coupling")
        subprocess.call(
            self.ppstr("postop_coupling_coeffs_fourrier_series"), shell=True
        )
        filename = self.tmp_dir + "/coupling_coeffs_fs.txt"
        tmp = femio.load_timetable(filename)
        return tmp.reshape((2 * self.M_fs + 1, self.neig))

    def get_spectral_elements(self, sort=False):
        eigval = self.postpro_eigenvalues()
        eigvect = self.postpro_eigenvectors()
        if sort:
            isort = np.argsort(eigval)
            eigval = eigval[isort]
            eigvect = eigvect[:, :, (isort)]
        return eigval, eigvect

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

    def postpro_modal_coupling(self):
        self.print_progress("QNMs coupling coeffs")
        subprocess.call(self.ppstr("postop_mode_coupling"), shell=True)
        filename = self.tmp_dir + "/mode_coupling.txt"
        tmp = femio.load_timetable(filename)
        return tmp  # .reshape((self.Ni_theta, self.neig))

    def postpro_modal_coupling_int(self):
        self.print_progress("QNMs coupling coeffs integrand")
        subprocess.call(self.ppstr("postop_mode_coupling_int"), shell=True)
        mode = femio.load_timetable(self.tmp_dir + "/mode_coupling_int.txt")
        u1 = np.zeros((self.Nix, self.Niy, self.neig), dtype=complex)
        u = mode.reshape((self.Niy, self.Nix, self.neig))
        for imode in range(self.neig):
            u1[:, :, imode] = np.flipud(u[:, :, imode]).T
        return u1

    def get_adjoint(self):
        return self.get_qty("adjoint.txt")

    def get_deq_deps(self):
        if self.pola is "TE":
            return self.get_qty("dEq_deps.txt")
        else:
            return self.get_qty("dEq_deps_x.txt"), self.get_qty("dEq_deps_y.txt")
