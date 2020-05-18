# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT

"""
Finite Element model of 2D media

"""

import numpy as np
from ..tools import femio
from ..basefem import *


class Scatt2D(BaseFEM):
    """A class for a finite element model of a 2D medium"""

    def __init__(self):
        super().__init__()
        self.dir_path = get_file_path(__file__)

        #: str: analysys type (either "direct" or "modal")
        self.analysis = "direct"
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
        self.theta_deg = 0.0

        #: line source position
        self.ls_flag = False
        self.xs = 0
        self.ys = 0

        #: beam?
        self.beam_flag = False
        self.waist = 1.5

        # opto-geometric parameters  -------------------------------------------
        self.h_pml = 2.0  #: flt: thickness pml
        self.hx_des = 2.0  #: flt: x - thickness scattering box (design)
        self.hy_des = 2.0  #: flt: y - thickness scattering box
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
        #: int: number of x interpolation points for near to far field calculations
        self.Nin2f_x = 500
        #: int: number of y interpolation points for near to far field calculations
        self.Nin2f_y = 500
        self.rat_n2f = 0.95
        self.inclusion_filename_ = "inclusion0.geo"

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

        self.yplane = 0

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
        return np.pi / 180.0 * (self.theta_deg)

    @property
    def omega0(self):
        return 2.0 * np.pi * self.cel / self.lambda0

    def _make_param_dict(self):
        param_dict = super()._make_param_dict()
        param_dict["ls_flag"] = int(self.ls_flag)
        param_dict["beam_flag"] = int(self.beam_flag)
        param_dict["target_flag"] = int(self.target_flag)
        return param_dict

    def get_field_ref_plane(self):
        self.postprocess("postop_fields_plane")
        u_tot = femio.load_table(self.tmppath("u_tot_ref_plane.out"))
        # u_i = femio.load_table(self.tmppath("u_i_point.txt"))
        # u = femio.load_table(self.tmppath("u_point.txt"))
        return u_tot

    def postpro_fields_box(self):
        path_T = self.tmppath("field_box_T.out")
        path_L = self.tmppath("field_box_L.out")
        path_R = self.tmppath("field_box_R.out")
        path_B = self.tmppath("field_box_B.out")
        self.postprocess("postop_fields_box")
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
        self.postprocess("postop_fields_n2f")
        u_out, vx_out, vy_out = {}, {}, {}
        for i in ["T", "L", "R", "B"]:
            u = femio.load_timetable(self.tmppath("field_n2f_{0}.out".format(i)))
            vx = femio.load_timetable(
                self.tmppath("field_dual_x_n2f_{0}.out".format(i))
            )
            vy = femio.load_timetable(
                self.tmppath("field_dual_y_n2f_{0}.out".format(i))
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
        k0 = 2 * np.pi / self.lambda0

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
        return np.abs(Itheta) ** 2 / (8 * np.pi)

    def get_field_map(self, name):
        field = femio.load_table(self.tmppath(name))
        return np.flipud(field.reshape((self.Niy, self.Nix)))

    def get_field_point(self):
        self.postprocess("postop_field_on_point")
        u_tot = femio.load_table(self.tmppath("u_tot_point.txt"))
        u_i = femio.load_table(self.tmppath("u_i_point.txt"))
        u = femio.load_table(self.tmppath("u_point.txt"))
        return u, u_tot, u_i

    def postpro_coupling_angle(self):
        self._print_progress("Angular sweep for coupling coeffs")
        self.postprocess("postop_coupling_coeffs_angle")
        filename = self.tmppath("coupling_coeffs.txt")
        tmp = femio.load_timetable(filename)
        return tmp.reshape((self.Ni_theta, self.neig))

    def postpro_fourrier_coefs_angle(self):
        self._print_progress("Fourrier coefficients for coupling")
        self.postprocess("postop_coupling_coeffs_fourrier_series")
        filename = self.tmppath("coupling_coeffs_fs.txt")
        tmp = femio.load_timetable(filename)
        return tmp.reshape((2 * self.M_fs + 1, self.neig))

    def postpro_modal_coupling(self):
        self._print_progress("QNMs coupling coeffs")
        self.postprocess("postop_mode_coupling")
        filename = self.tmppath("mode_coupling.txt")
        tmp = femio.load_timetable(filename)
        return tmp  # .reshape((self.Ni_theta, self.neig))

    def postpro_modal_coupling_int(self):
        self._print_progress("QNMs coupling coeffs integrand")
        self.postprocess("postop_mode_coupling_int")
        mode = femio.load_timetable(self.tmppath("mode_coupling_int.txt"))
        u1 = np.zeros((self.Nix, self.Niy, self.neig), dtype=complex)
        u = mode.reshape((self.Niy, self.Nix, self.neig))
        for imode in range(self.neig):
            u1[:, :, imode] = np.flipud(u[:, :, imode]).T
        return u1

    def get_deq_deps(self):
        if self.pola is "TE":
            return self._get_qty("dEq_deps.txt")
        else:
            return self._get_qty("dEq_deps_x.txt"), self._get_qty("dEq_deps_y.txt")
