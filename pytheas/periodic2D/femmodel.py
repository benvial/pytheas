# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


import os
import numpy as np
import scipy as sc
from ..tools import femio
from ..basefem import *


class Periodic2D(BaseFEM):
    """A class for a finite element model of a 2D mono-periodic medium.

    The model consist of a single unit cell with quasi-periodic boundary conditions
    in the :math:`x` direction enclosed with perfectly matched layers (PMLs)
    in the :math:`y` direction to truncate the semi infinite media. From top to bottom:

    - PML top
    - superstrate (incident medium)
    - layer 2
    - design layer: this is the layer containing the periodic pattern, can be continuous or discrete
    - layer 1
    - substrate
    - PML bottom

    Parameters
    ----------
    analysis : str, default "direct"
        Analysis type: either "direct" (plane wave) or
        "modal" (spectral problem)

    pola : str, default "TE"
        Polarization case: either "TE" (E along z) or "TM" (H along z)

    A : float, default 1
        Incident plane wave amplitude

    lambda0 : float, default 1
        Incident plane wave wavelength in free space

    lambda_mesh : float, default 1
        Wavelength to use for meshing

    theta_deg : float, default 0
        Incident plane wave angle (in degrees).
        Light comes from the top (travels along -y if normal incidence,
        ``theta_deg=0`` is set)

    d : float, default 0.8
        Periodicity

    h_sup : float, default 1
        Thickness superstrate

    h_sub : float, default 1
        Thickness substrate

    h_layer1 : float, default 0.1
        Thickness layer 1

    h_layer2 : float, default 0.1
        Thickness layer 2

    h_des : float, default 1
        Thickness layer design

    h_pmltop : float, default 1
        Thickness pml top

    h_pmlbot : float, default 1
        Thickness pml bot

    a_pml : float, default 1
        PMLs complex y-stretching parameter, real part

    b_pml : float, default 1
        PMLs complex y-stretching parameter, imaginary part

    eps_sup : complex, default (1 - 0 * 1j)
        Permittivity superstrate

    eps_sub : complex, default (1 - 0 * 1j)
        Permittivity substrate

    eps_layer1 : complex, default (1 - 0 * 1j)
        Permittivity layer 1

    eps_layer2 : complex, default (1 - 0 * 1j)
        Permittivity layer 2

    eps_des : complex, default (1 - 0 * 1j)
        Permittivity layer design

    eps_incl : complex, default (1 - 0 * 1j)
        Permittivity inclusion

    """

    def __init__(
        self,
        analysis="direct",
        pola="TE",
        A=1,
        lambda0=1,
        lambda_mesh=1,
        theta_deg=0,
        d=0.8,
        h_sup=1,
        h_sub=1,
        h_layer1=0.1,
        h_layer2=0.1,
        h_des=1.0,
        h_pmltop=1.0,
        h_pmlbot=1.0,
        a_pml=1,
        b_pml=1,
        eps_sup=1 - 0 * 1j,
        eps_sub=1 - 0 * 1j,
        eps_layer1=1 - 0 * 1j,
        eps_layer2=1 - 0 * 1j,
        eps_des=1 - 0 * 1j,
        eps_incl=1 - 0 * 1j,
        mu_incl=1 - 0 * 1j,
        mu_des=1 - 0 * 1j,
    ):

        super().__init__()
        self.dir_path = get_file_path(__file__)
        self.analysis = analysis
        self.pola = pola
        self.A = A
        self.lambda0 = lambda0
        self.lambda_mesh = lambda_mesh
        self.theta_deg = theta_deg
        # opto-geometric parameters  -------------------------------------------
        self.d = d
        self.h_sup = h_sup
        self.h_sub = h_sub
        self.h_layer1 = h_layer1
        self.h_layer2 = h_layer2
        self.h_des = h_des
        self.h_pmltop = h_pmltop
        self.h_pmlbot = h_pmlbot
        self.a_pml = a_pml
        self.b_pml = b_pml
        self.eps_sup = eps_sup
        self.eps_sub = eps_sub
        self.eps_layer1 = eps_layer1
        self.eps_layer2 = eps_layer2
        self.eps_des = eps_des
        self.eps_incl = eps_incl
        self.mu_incl = mu_incl
        self.mu_des = mu_des
        self.dom_des = 4000
        self.nb_incl = 1
        self.aniso = False
        # if not isinstance(self.eps_layer1, complex):
        #     raise TypeError("bar must be a complex number")

        # postprocessing -------------------------------------------------
        # int: number of diffraction orders
        # for postprocessing diffraction efficiencies
        # N_d_order = 0
        # int: number of x integration points
        # for postprocessing diffraction efficiencies
        self.npt_integ = 1000
        # int: number of y slices points
        # for postprocessing diffraction efficiencies
        self.nb_slice = 20
        # flt:  such that `scan_dist  = min(h_sup, h_sub)/scan_dist_ratio`
        self.scan_dist_ratio = 5
        self.delta_x, self.delta_y = 0, 0

        self.nper = 3
        # modal analysis -------------------------------------------------
        # int: number of eigenvalues searched for in modal analysis
        self.neig = 6
        # flt: wavelength around which to search eigenvalues
        self.lambda0search = 1

        # topology optimization   ------------------------
        self.adjoint = False
        self.r_tar = 1 + 0j
        self.t_tar = 0 + 0j

    @property
    def scan_dist(self):
        return min(self.h_sub, self.h_sup) / self.scan_dist_ratio

    @property
    def ycut_sub_min(self):
        return -self.h_sub + self.scan_dist

    @property
    def ycut_sub_max(self):
        return -self.scan_dist

    @property
    def ycut_sup_min(self):
        return self.h_layer1 + self.h_layer2 + self.h_des + self.scan_dist

    @property
    def ycut_sup_max(self):
        return self.h_layer1 + self.h_layer2 + self.h_des + self.h_sup - self.scan_dist

    @property
    def domX_L(self):
        return -self.d / 2

    @property
    def domX_R(self):
        return self.d / 2

    @property
    def domY_B(self):
        return -self.h_sub + self.delta_y  # - self.h_pmlbot

    @property
    def domY_T(self):
        # + self.h_pmltop
        return self.h_layer1 + self.h_layer2 + self.h_des + self.h_sup - self.delta_y

    @property
    def corners_des(self):
        return -self.d / 2, +self.d / 2, self.h_layer1, self.h_layer1 + self.h_des

    @property
    def theta(self):
        return np.pi / 180.0 * (self.theta_deg)

    @property
    def omega0(self):
        return 2.0 * np.pi * self.cel / self.lambda0

    @property
    def N_d_order(self):
        a = self.d / self.lambda0
        b = np.sin(self.theta)
        WRA = np.zeros((4))
        WRA[0] = a * (np.sqrt(self.eps_sup.real) + b)
        WRA[1] = a * (np.sqrt(self.eps_sub.real) + b)
        WRA[2] = a * (np.sqrt(self.eps_sup.real) - b)
        WRA[3] = a * (np.sqrt(self.eps_sub.real) - b)
        return int(max(abs(WRA)))

    # def _make_param_dict(self):
    #     param_dict = super()._make_param_dict()
    #     param_dict["aniso"] = int(self.aniso)
    #     return param_dict

    def compute_solution(self):
        """Compute the solution of the FEM problem using getdp"""
        res_list = ["helmoltz_scalar", "helmoltz_scalar_modal"]
        super().compute_solution(res_list=res_list)

    def postpro_absorption(self):
        """ Compute the absorption coefficient

        Returns
        -------
        Q : float
            Absorption coefficient
        """
        self._print_progress("Postprocessing absorption")
        self.postprocess("postop_absorption")
        Q = femio.load_table(self.tmppath("Q.txt")).real
        return Q

    def _postpro_fields_cuts(self):
        """ Compute the field cuts in substrate and superstarte

        Returns
        -------
        u_diff_t : array-like
            Transmitted field cuts
        u_diff_r : array-like
            Reflected field cuts
        """
        path_sub = self.tmppath("sub_field_cuts.out")
        path_sup = self.tmppath("sup_field_cuts.out")
        if os.path.isfile(path_sub):
            os.remove(path_sub)
        if os.path.isfile(path_sup):
            os.remove(path_sup)
        self.postprocess("postop_fields_cuts")
        u_diff_t = femio.load_table(path_sub)
        u_diff_r = femio.load_table(path_sup)
        return u_diff_t, u_diff_r

    def get_field_map(self, name):
        """Retrieve a field map.

        Parameters
        ----------
        name : str
            Choose between "u" (scattered field), "u_tot" (total field)

        Returns
        -------
        array
            A 2D complex array of shape (`Nix`, `Niy`)

        """
        field = femio.load_table(self.tmppath(name + ".txt"))
        field = np.flipud(field.reshape((self.Niy, self.Nix))).T
        return field

    def diffraction_efficiencies(self, cplx_effs=False, orders=False):
        """Postprocess diffraction efficiencies.

        Parameters
        ----------
        cplx_effs : bool
            If `True`, return complex coefficients (amplitude reflection and transmission).
            If `False`, return real coefficients (power reflection and transmission)
        orders : bool
            If `True`, computes the transmission and reflection for all the propagating diffraction orders.
            If `False`, returns the sum of all the propagating diffraction orders.

        Returns
        -------
        dict
            A dictionary containing the diffraction efficiencies.

        """
        self._print_progress("Processing diffraction efficiencies")
        if self.pola == "TE":
            nu = 1
        else:
            nu = 1 / self.eps_sub
        order_shift = 0
        No_ordre = np.linspace(
            -self.N_d_order + order_shift,
            self.N_d_order + order_shift,
            2 * self.N_d_order + 1,
        )
        x_slice = np.linspace(-self.d / 2, self.d / 2, self.npt_integ)
        y_slice_t = np.linspace(self.ycut_sub_min, self.ycut_sub_max, self.nb_slice)
        y_slice_r = np.linspace(self.ycut_sup_min, self.ycut_sup_max, self.nb_slice)
        k_sub = 2 * np.pi * np.sqrt(self.eps_sub) / self.lambda0
        k_sup = 2 * np.pi * np.sqrt(self.eps_sup) / self.lambda0
        alpha_sup = -k_sup * np.sin(self.theta)
        beta_sup = np.sqrt(k_sup ** 2 - alpha_sup ** 2)
        # beta_sub = np.sqrt(k_sub ** 2 - alpha_sup ** 2)
        s_t = np.zeros((1, (2 * self.N_d_order + 1)), complex)[0, :]
        s_r = np.zeros((1, (2 * self.N_d_order + 1)), complex)[0, :]
        Aeff_t = np.zeros((self.nb_slice, 2 * self.N_d_order + 1), complex)
        Aeff_r = np.zeros((self.nb_slice, 2 * self.N_d_order + 1), complex)

        field_diff_t, field_diff_r = self._postpro_fields_cuts()

        field_diff_t = np.transpose(
            field_diff_t.reshape(self.npt_integ, self.nb_slice, order="F")
        )
        field_diff_r = np.transpose(
            field_diff_r.reshape(self.npt_integ, self.nb_slice, order="F")
        )
        alphat_t = alpha_sup + 2 * np.pi / (self.d) * No_ordre
        alphat_r = alpha_sup + 2 * np.pi / (self.d) * No_ordre
        betat_sup = np.conj(np.sqrt(k_sup ** 2 - alphat_r ** 2))
        betat_sub = np.conj(np.sqrt(k_sub ** 2 - alphat_t ** 2))
        for m1 in range(0, self.nb_slice):
            slice_t = field_diff_t[m1, :]
            slice_r = field_diff_r[m1, :]

            for k in range(0, 2 * self.N_d_order + 1):
                expalpha_t = np.exp(-1j * alphat_t[k] * x_slice)
                expalpha_r = np.exp(-1j * alphat_r[k] * x_slice)
                s_t[k] = np.trapz(slice_t * expalpha_t, x=x_slice) / self.d
                s_r[k] = np.trapz(slice_r * expalpha_r, x=x_slice) / self.d

            Aeff_t[m1, :] = s_t * np.exp(-1j * betat_sub * (y_slice_t[m1]))
            Aeff_r[m1, :] = s_r * np.exp(
                1j * betat_sup * (y_slice_r[m1] - (self.ycut_sup_min - self.scan_dist))
            )

        # Aeff_r = -np.conj(Aeff_r)

        Beff_t = (np.abs(Aeff_t)) ** 2 * betat_sub / beta_sup * nu
        Beff_r = (np.abs(Aeff_r)) ** 2 * betat_sup / beta_sup

        # print(Aeff_r)
        # print(Aeff_t)

        rcplx = np.mean(Aeff_r, axis=0)
        tcplx = np.mean(Aeff_t, axis=0)

        Rorders = np.mean(Beff_r.real, axis=0)
        Torders = np.mean(Beff_t.real, axis=0)
        R = np.sum(Rorders, axis=0)
        T = np.sum(Torders, axis=0)
        Q = self.postpro_absorption()
        B = T + R + Q  # energy balance
        if self.python_verbose:
            print("  Energy balance")
            print(
                "    R            = ",
                "%0.6f" % R,
                "     (std slice2slice =",
                "%0.6e" % np.std(np.sum(Aeff_r.real, axis=1)),
                ")",
            )
            print(
                "    T            = ",
                "%0.6f" % T,
                "     (std slice2slice =",
                "%0.6e" % np.std(np.sum(Aeff_t.real, axis=1)),
                ")",
            )
            print("    Q            = ", "%0.6f" % Q)
            print("    ------------------------")
            print("    Total        = ", "%0.6f" % (B))

        if cplx_effs:
            R, T = rcplx, tcplx
        else:
            if orders:
                R, T = Rorders, Torders
        eff = dict()
        eff["R"] = R
        eff["T"] = T
        eff["Q"] = Q
        eff["B"] = B
        return eff

    def pseudo_periodize(self, qt):
        k0 = 2 * np.pi / self.lambda0
        coef = np.exp(-1j * k0 * np.sin(self.theta) * self.d)
        qtper = qt
        for i in range(self.nper - 1):
            qtper = np.row_stack((qtper, qt * coef ** (i + 1)))
        return qtper

    def periodize(self, qt):
        qtper = qt
        for _ in range(self.nper - 1):
            qtper = np.row_stack((qtper, qt))
        return qtper

    def plot_field_and_pattern(
        self,
        fig,
        ax,
        field,
        pattern,
        cmap_div,
        cmap_mat,
        cbar=True,
        vmin=None,
        vmax=None,
    ):

        self._print_progress("Plotting field map")
        # x = np.linspace(
        #     self.nper * self.domX_L, self.nper * self.domX_R, self.nper * self.Nix
        # )
        # y = np.linspace(self.domY_B, self.domY_T, self.Niy)

        # xx, yy = np.meshgrid(x, y)
        extent = (
            self.nper * self.domX_L,
            self.nper * self.domX_R,
            self.domY_B,
            self.domY_T,
        )

        im1 = ax.imshow(
            field.T,
            interpolation="bilinear",
            cmap=cmap_div,
            vmin=vmin,
            vmax=vmax,
            extent=extent,
        )
        if cbar:
            fig.colorbar(im1, fraction=0.046, pad=0.04)
        ax.imshow(
            pattern.T,
            interpolation="None",
            cmap=cmap_mat,
            alpha=0.4,
            extent=(
                self.nper * self.domX_L,
                self.nper * self.domX_R,
                self.h_layer1,
                self.h_layer1 + self.h_des,
            ),
        )
        ax.imshow(
            field.T,
            alpha=0.0,
            interpolation="bilinear",
            cmap=cmap_div,
            vmin=vmin,
            vmax=vmax,
            extent=(
                self.nper * self.domX_L,
                self.nper * self.domX_R,
                self.domY_B,
                self.domY_T,
            ),
        )
        ax.set_ylim((self.domY_B, self.domY_T))

    def plot_pattern_per(self, fig, ax, pattern, nper, cmap):
        im1 = ax.imshow(
            pattern.T,
            cmap=cmap,
            extent=(
                -self.nper * self.d / 2,
                self.nper * self.d / 2,
                self.h_layer1,
                self.h_layer1 + self.h_des,
            ),
        )
        fig.colorbar(im1, fraction=0.046, pad=0.04)

    def plot_fieldv_streamlines(self, ax, vx, vy, cmap):
        self._print_progress("Plotting vector field streamlines")
        x = np.linspace(
            self.nper * self.domX_L, self.nper * self.domX_R, self.nper * self.Nix
        )
        y = np.linspace(self.domY_B, self.domY_T, self.Niy)
        vlognorm = (np.sqrt(np.abs(vx) ** 2 + np.abs(vy) ** 2)).T
        Ndens = 1
        density = [Ndens * self.nper, Ndens * self.Niy / self.Nix]
        lw = 3 * vlognorm / vlognorm.max()
        ax.streamplot(
            x,
            y,
            vx.T,
            vy.T,
            color="k",
            linewidth=lw,
            cmap=cmap,
            density=density,
            arrowstyle="Fancy",
        )
        ax.set_ylim((self.domY_B, self.domY_T))
