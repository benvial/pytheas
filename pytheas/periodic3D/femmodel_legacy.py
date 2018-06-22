import numpy as np
import scipy as sc
import os
import subprocess
import shutil
import pytheas.tools.femio as femio
# dir_path = os.path.dirname(os.path.abspath(__file__))
pi = np.pi


class PeriodicMediumFEM3D:
    """A class for a finite element model of a 3D bi-periodic
       medium using Gmsh_ and GetDP_.

        .. _Gmsh:
            http://gmsh.info/
        .. _GetDP:
            http://getdp.info/
    """

    dir_path = os.path.dirname(os.path.abspath(__file__))
    ID = "fem"  #: str: name of simulation
    #: str: analisys type (either "diffraction" or "modal")
    analysis = "diffraction"
    # filenames
    geom_filename = os.path.join(
        dir_path, "base/geometry.geo")  #: str: Gmsh geometry filename
    pro_filename = os.path.join(dir_path,
                                "base/main.pro")  #: str: GetDP pro filename
    content_geo = femio.get_content(geom_filename)
    content_pro = femio.get_content(pro_filename)

    tmp_dir = "./tmp"
    path_pos = None

    el_order = 1
    celltype = "tetra"
    type_des = "nodes"

    getdp_verbose = 0  #: str: GetDP verbose (int between 0 and 4)
    gmsh_verbose = 0  #: str: Gmsh verbose (int between 0 and 4)
    python_verbose = 1  #: str: python verbose (int between 0 and 1)
    epsilon0 = 8.854187817e-12  #: flt: vacuum permittivity
    mu0 = 4. * pi * 1e-7  #: flt: vacuum permeability
    #: flt: speed of light in vacuum
    cel = 1.0 / (np.sqrt(epsilon0 * mu0))

    # Incident plane wave parameters
    #: flt: incident plane wave amplitude
    A = 1.0
    #: flt: incident plane wave wavelength in free space
    lambda0 = 1.0
    #: flt : incident plane wave angle (in degrees).
    #: Light comes from the top
    #: (travels along -y if normal incidence, `theta_deg=0` is set)
    theta_deg = 0.
    phi_deg = 0
    psi_deg = 0
    Deph = 0

    #: flt: global mesh parameter
    #: `MeshElementSize = lambda0/(parmesh*n)`, `n`: refractive index
    parmesh = 10.
    #: flt: design subdomain mesh parameter
    parmesh_des = 10.
    #: flt: PMLs mesh parameter
    parmesh_pml = 7.

    # opto-geometric parameters  -------------------------------------------
    #: flt: periods
    period_x = 1
    period_y = 1
    thick_L1 = 0.1  #: flt: thickness layer 1 (superstrate)
    thick_L2 = 0.1  #: flt: thickness layer 2
    thick_L3 = 0.1  #: flt: thickness layer 3 (interp)
    thick_L4 = 0.1  #: flt: thickSness layer 4
    thick_L5 = 0.1  #: flt: thickness layer 5
    thick_L6 = 0.1  #: flt: thickness layer 6 (substrate)
    PML_top = 1.  #: flt: thickness pml top
    PML_bot = 1.  #: flt: thickness pml bot
    #: flt: PMLs parameter, real part
    a_pml = 1  #: flt: PMLs parameter, real part
    b_pml = 1  #: flt: PMLs parameter, imaginary part
    eps_L1 = 1  #: flt: permittivity layer 1 (superstrate)
    eps_L2 = 1  #: flt: permittivity layer 2
    eps_L3 = 1  #: flt: permittivity layer 3
    eps_L4 = 1  #: flt: permittivity layer 4
    eps_L5 = 1  #: flt: permittivity layer 5
    eps_L6 = 1  #: flt: permittivity layer 6 (substrate)

    dom_des = 5000  #: design domain number (check .geo/.pro files)

    # postprocessing -------------------------------------------------
    #: int: number of diffraction orders
    #: for postprocessing diffraction efficiencies
    N_d_order = 0
    orders = False
    cplx_effs = False
    eff_verbose = False
    #: int: number of x integration points
    #: for postprocessing diffraction efficiencies
    ninterv_integ = 60
    #: int: number of y slices points
    #: for postprocessing diffraction efficiencies
    nb_slice = 5
    #: flt:  such that `scan_dist  = min(h_sup, hsub)/scan_dist_ratio`
    scan_dist_ratio = 5

    @property
    def zmin_interp(self):
        return self.thick_L5 + self.thick_L4

    @property
    def zmax_interp(self):
        return self.zmin_interp + self.thick_L3

    @property
    def scan_dist(self):
        return min(self.thick_L1, self.thick_L6) / self.scan_dist_ratio

    @property
    def theta_0(self):
        return pi / 180. * (self.theta_deg)

    @property
    def phi_0(self):
        return pi / 180. * (self.phi_deg)

    @property
    def psi_0(self):
        return pi / 180. * (self.psi_deg)

    @property
    def omega0(self):
        return 2. * pi * self.cel / self.lambda0

    @property
    def param_dict(self):
        return self.make_param_dict()

    @property
    def path_geo(self):
        return os.path.join(self.tmp_dir, "geometry.geo")

    @property
    def path_pro(self):
        return os.path.join(self.tmp_dir, "main.pro")

    @property
    def path_mesh(self):
        return os.path.join(self.tmp_dir, "mesh.msh")

    @property
    def content_par(self):
        return femio.make_inputs(self.param_dict)

    # @property
    # def N_d_order(self):
    #     N = self.d/self.lambda0 * (np.sqrt([self.eps_L1, self.eps_L6]) - np.sin(self.theta))
    #     return int(max(N))

    def initialize(self):
        try:
            os.mkdir(self.tmp_dir)
        except FileExistsError:
            pass
            # shutil.rmtree(tmp_dir)
        # create tmp parameters files files
        femio.maketmp(self.content_par, "parameters.dat", dirname=self.tmp_dir)
        # create tmp geo file
        femio.maketmp(self.content_geo, "geometry.geo", dirname=self.tmp_dir)
        # create tmp geo file
        femio.maketmp(self.content_pro, "main.pro", dirname=self.tmp_dir)

    def update_params(self):
        femio.maketmp(self.content_par, "parameters.dat", dirname=self.tmp_dir)

    def cleanup(self):
        os.remove("rm *.msh *.pre *.res *.dat *.txt *.pyc")

    #
    # def make_tmp_files(self):
    #     path_pro = femio.maketmp(pro_merged, suffix=".pro", dir=self.tmp_dir)

    def make_param_dict(self):
        layer_diopter = self.ancillary_problem()
        param_dict = dict()
        param_dict["parmesh"] = self.parmesh
        param_dict["parmesh_des"] = self.parmesh_des
        param_dict["parmesh_pml"] = self.parmesh_pml
        param_dict["el_order"] = self.el_order
        param_dict["lambda0"] = self.lambda0
        param_dict["theta0"] = self.theta_0
        param_dict["phi0"] = self.phi_0
        param_dict["psi0"] = self.psi_0
        param_dict["cel"] = self.cel
        param_dict["mu0"] = self.mu0
        param_dict["epsilon0"] = self.epsilon0
        param_dict["a_pml"] = self.a_pml
        param_dict["b_pml"] = self.b_pml
        param_dict["period_x"] = self.period_x
        param_dict["period_y"] = self.period_y
        nb_layer = 6
        layer = []
        for k1 in range(0, nb_layer):
            layer.append({})
        layer[0]['epsilon'] = self.eps_L1
        layer[1]['epsilon'] = self.eps_L2
        layer[2]['epsilon'] = self.eps_L3
        layer[3]['epsilon'] = self.eps_L4
        layer[4]['epsilon'] = self.eps_L5
        layer[5]['epsilon'] = self.eps_L6

        layer[0]['thickness'] = self.thick_L1
        layer[1]['thickness'] = self.thick_L2
        layer[2]['thickness'] = self.thick_L3
        layer[3]['thickness'] = self.thick_L4
        layer[4]['thickness'] = self.thick_L5
        layer[5]['thickness'] = self.thick_L6

        layer[nb_layer - 2]['hh'] = 0
        layer[nb_layer
              - 1]['hh'] = layer[nb_layer - 2]['hh'] - layer[nb_layer
                                                             - 1]['thickness']
        for k in range(nb_layer - 3, -1, -1):
            layer[k]['hh'] = layer[k + 1]['hh'] + layer[k + 1]['thickness']
        for i5 in range(0, nb_layer):
            param_dict['eps_re_L_' + str(i5 + 1)] = layer[i5]['epsilon'].real
            param_dict['eps_im_L_' + str(i5 + 1)] = layer[i5]['epsilon'].imag
            param_dict['thick_L_' + str(i5 + 1)] = layer[i5]['thickness']
            param_dict['hh_L_' + str(i5 + 1)] = layer[i5]['hh']
        param_dict['PML_bot'] = self.PML_bot
        param_dict['PML_top'] = self.PML_top
        param_dict['PML_bot_hh'] = layer[-1]['hh'] - self.PML_bot
        param_dict['PML_top_hh'] = layer[0]['hh'] + self.thick_L1
        param_dict['Expj_subs_re'] = layer_diopter[1]['Psi'][0].real
        param_dict['Exmj_subs_re'] = layer_diopter[1]['Psi'][1].real
        param_dict['Eypj_subs_re'] = layer_diopter[1]['Psi'][2].real
        param_dict['Eymj_subs_re'] = layer_diopter[1]['Psi'][3].real
        param_dict['Ezpj_subs_re'] = layer_diopter[1]['Psi'][4].real
        param_dict['Ezmj_subs_re'] = layer_diopter[1]['Psi'][5].real
        param_dict['Expj_subs_im'] = layer_diopter[1]['Psi'][0].imag
        param_dict['Exmj_subs_im'] = layer_diopter[1]['Psi'][1].imag
        param_dict['Eypj_subs_im'] = layer_diopter[1]['Psi'][2].imag
        param_dict['Eymj_subs_im'] = layer_diopter[1]['Psi'][3].imag
        param_dict['Ezpj_subs_im'] = layer_diopter[1]['Psi'][4].imag
        param_dict['Ezmj_subs_im'] = layer_diopter[1]['Psi'][5].imag
        param_dict['gamma_subs_re'] = layer_diopter[1]['gamma'].real
        param_dict['gamma_subs_im'] = layer_diopter[1]['gamma'].imag

        param_dict['Expj_super_re '] = layer_diopter[0]['Psi'][0].real
        param_dict['Exmj_super_re '] = layer_diopter[0]['Psi'][1].real
        param_dict['Eypj_super_re '] = layer_diopter[0]['Psi'][2].real
        param_dict['Eymj_super_re '] = layer_diopter[0]['Psi'][3].real
        param_dict['Ezpj_super_re '] = layer_diopter[0]['Psi'][4].real
        param_dict['Ezmj_super_re '] = layer_diopter[0]['Psi'][5].real
        param_dict['Expj_super_im '] = layer_diopter[0]['Psi'][0].imag
        param_dict['Exmj_super_im '] = layer_diopter[0]['Psi'][1].imag
        param_dict['Eypj_super_im '] = layer_diopter[0]['Psi'][2].imag
        param_dict['Eymj_super_im '] = layer_diopter[0]['Psi'][3].imag
        param_dict['Ezpj_super_im '] = layer_diopter[0]['Psi'][4].imag
        param_dict['Ezmj_super_im '] = layer_diopter[0]['Psi'][5].imag
        param_dict['gamma_super_re '] = layer_diopter[0]['gamma'].real
        param_dict['gamma_super_im '] = layer_diopter[0]['gamma'].imag

        param_dict["ninterv_integ"] = self.ninterv_integ
        param_dict["nb_slice"] = self.nb_slice
        param_dict["scan_dist"] = self.scan_dist

        return param_dict

    def get_design_nodes(self):
        return femio.get_nodes(self.path_mesh, self.dom_des, self.celltype)

    def get_design_elements(self):
        return femio.get_elements(self.path_mesh, self.dom_des, self.celltype)

    def make_eps_pos(self, des_ID, eps_des):
        # create a pos file to be read by getdp
        eps_des_pos = femio.make_pos(des_ID, eps_des, self.content_mesh,
                                     "eps_des", type=self.type_des)
        return femio.maketmp(eps_des_pos, "eps_des.pos", dirname=self.tmp_dir)

    def make_mesh_pos(self, els, nodes):
        return femio.make_content_mesh_pos(nodes, els, self.dom_des, self.celltype)

    def make_mesh(self, other_option=""):
        femio.mesh_model(
            self.path_mesh,
            self.path_geo,
            dim=[1, 2, 3],
            verbose=self.gmsh_verbose,
            other_option=other_option)
        content_mesh = femio.get_content(self.path_mesh)
        return content_mesh

    def compute_solution(self, **kwargs):
        if self.analysis == "diffraction":
            argstr = "-petsc_prealloc 500 -ksp_type preonly \
                     -pc_type lu -pc_factor_mat_solver_package mumps"

            resolution = "helmholtz_vector"
        elif self.analysis == "modal":
            argstr = "-slepc -eps_type krylovschur \
                       -st_ksp_type preonly \
                       -st_pc_type lu \
                       -st_pc_factor_mat_solver_package mumps \
                       -eps_max_it 300 \
                       -eps_target 0.00001 \
                       -eps_target_real \
                       -eps_mpd 600 -eps_nev 400"

            resolution = "helmholtz_vector_modal"
        else:
            raise TypeError(
                "Wrong analysis specified: choose between diffraction and modal"
            )

        femio.solve_problem(
            resolution,
            self.path_pro,
            self.path_mesh,
            verbose=self.getdp_verbose,
            path_pos=self.path_pos,
            argstr=argstr)

    def make_fdens(self, pattern):
        n_x, n_y, n_z = pattern.shape
        x = np.linspace(-self.period_x / 2, self.period_x / 2, n_x + 1)
        y = np.linspace(-self.period_y / 2, self.period_y / 2, n_y + 1)
        z = np.linspace(self.zmin_interp, self.zmax_interp, n_z + 1)
        dx, dy, dz = x[1] - x[0], y[1] - y[0], z[1] - z[0]
        x0, x1 = x[0], x[-1]
        y0, y1 = y[0], y[-1]
        z0, z1 = z[0], z[-1]
        x = np.linspace(x0 + dx / 2, x1 - dx / 2, n_x)
        y = np.linspace(y0 + dy / 2, y1 - dy / 2, n_y)
        z = np.linspace(z0 + dz / 2, z1 - dz / 2, n_z)
        xx, yy, zz = np.meshgrid(x, y, z)
        points = np.vstack((xx.ravel(), yy.ravel(), zz.ravel())).T
        pat = pattern.reshape(n_x * n_y * n_z)
        fdens = sc.interpolate.NearestNDInterpolator(points, pat)
        return fdens

    def assign_material(self, mat, matprop, density, lambda0):
        pattern = mat.mat_rand
        eps_nodes = np.zeros_like(density, dtype=complex)
        eps_pattern = np.zeros_like(pattern, dtype=complex)
        for i in range(mat.nb_threshold):
            if isinstance(matprop[i], str):
                ncomplex = ri.get_complex_index(lambda0, matprop[i])
            else:
                ncomplex = matprop[i]
            eps_nodes[density == mat.threshold_val[i]] = ncomplex**2
            eps_pattern[pattern == mat.threshold_val[i]] = ncomplex**2
        return eps_nodes, eps_pattern

    def ppstr(self, postop):
        return femio.postprostring(postop, self.path_pro, self.path_mesh,
                                   self.path_pos, self.getdp_verbose)

    def postpro_choice(self, name, filetype):
        if filetype in {"pos", "txt"}:
            subprocess.call(self.ppstr(name + "_" + filetype), shell=True)
        else:
            raise TypeError(
                "Wrong filetype specified: choose between txt and pos")

    def postpro_absorption(self):
        subprocess.call(self.ppstr("postopQ"), shell=True)
        path = self.tmp_dir + "/Q.txt"
        Q = np.loadtxt(
            path, skiprows=0, usecols=[1]) + 1j * np.loadtxt(
                path, skiprows=0, usecols=[1])
        return Q.real

    def postpro_fields_cuts(self):
        subprocess.call(self.ppstr("Ed") + " -order 2", shell=True)

    def postpro_epsilon(self):
        subprocess.call(
            [self.ppstr("postop_epsilon") + " -order 2"], shell=True)

    def postpro_fields(self, filetype="txt"):
        self.postpro_choice("postop_fields", filetype)

    def get_field_map(self, name):
        field = femio.load_table(self.tmp_dir + "/" + name)
        return field.reshape((self.Niy, self.Nix)).T

    def get_objective(self, run=False):
        if not self.adjoint:
            subprocess.call(self.ppstr("postop_int_objective"), shell=True)
        return femio.load_table(self.tmp_dir + "/objective.txt").real

    def get_adjoint(self):
        return femio.load_node_table(self.tmp_dir + "/adjoint.txt")

    def get_deq_deps(self):
        return femio.load_node_table(self.tmp_dir + "/dEq_deps.txt")

    def postpro_eigenvalues(self):
        subprocess.call(self.ppstr("postop_eigenvalues"), shell=True)
        filename = self.tmp_dir + "/EigenValues.txt"
        re = np.loadtxt(filename, usecols=[1])
        im = np.loadtxt(filename, usecols=[5])
        return re + 1j * im

    def postpro_eigenvectors(self, filetype="txt"):
        self.postpro_choice("postop_eigenvectors", filetype)
        if filetype is "txt":
            ev = femio.load_timetable(self.tmp_dir + "/EigenVectors.txt")
            return ev.reshape((self.Nix, self.Niy, self.neig))

    def diffraction_efficiencies(self):
        self.postpro_fields_cuts()
        npt_integ = self.ninterv_integ + 1
        # print('gmsh cuts done !')
        period_x, period_y = self.period_x, self.period_y
        N_d_order = self.N_d_order
        lambda0 = self.lambda0
        theta_0 = self.theta_0
        phi_0 = self.phi_0
        nb_slice = self.nb_slice
        x_t = np.linspace(-period_x / 2, period_x / 2, npt_integ)
        x_r = x_t
        y_t = np.linspace(-period_y / 2, period_y / 2, npt_integ)
        y_r = y_t
        decalage = 0
        No_ordre = np.linspace(-N_d_order + decalage, N_d_order + decalage,
                               2 * N_d_order + 1)
        Nb_ordre = No_ordre.shape[0]
        alpha0 = 2 * pi / lambda0 * np.sin(theta_0) * np.cos(phi_0)
        beta0 = 2 * pi / lambda0 * np.sin(theta_0) * np.sin(phi_0)
        gamma0 = 2 * pi / lambda0 * np.cos(theta_0)
        alphat = alpha0 + 2 * pi / period_x * No_ordre
        betat = beta0 + 2 * pi / period_y * No_ordre
        gammatt = np.zeros((Nb_ordre, Nb_ordre), dtype=complex)
        gammatr = np.zeros((Nb_ordre, Nb_ordre), dtype=complex)
        AXsir = np.zeros((Nb_ordre, Nb_ordre, nb_slice), dtype=complex)
        AXsit = np.zeros((Nb_ordre, Nb_ordre, nb_slice), dtype=complex)

        nb_layer_diopter = 2
        layer_diopter = []
        for k1 in range(0, nb_layer_diopter):
            layer_diopter.append({})
        layer_diopter[0]['epsilon'] = self.eps_L1
        layer_diopter[1]['epsilon'] = self.eps_L6
        layer_diopter[0]['kp'] = 2 * pi / lambda0 * sc.sqrt(
            layer_diopter[0]['epsilon'])
        layer_diopter[1]['kp'] = 2 * pi / lambda0 * sc.sqrt(
            layer_diopter[1]['epsilon'])
        layer_diopter[0]['gamma'] = sc.sqrt(layer_diopter[0]['kp']**2 - alpha0
                                            ** 2 - beta0**2)
        layer_diopter[1]['gamma'] = sc.sqrt(layer_diopter[1]['kp']**2 - alpha0
                                            ** 2 - beta0**2)

        for nt in range(0, Nb_ordre):
            for mt in range(0, Nb_ordre):
                gammatt[nt, mt] = np.sqrt(layer_diopter[-1]['kp']**2 -
                                          alphat[nt]**2 - betat[mt]**2)
        for nr in range(0, Nb_ordre):
            for mr in range(0, Nb_ordre):
                gammatr[nr, mr] = np.sqrt(layer_diopter[0]['kp']**2 -
                                          alphat[nr]**2 - betat[mr]**2)

        path_t = self.tmp_dir + "/Etot_XYcut.out"
        path_r = self.tmp_dir + "/Edif_XYcut.out"
        Ex_t2 = np.loadtxt(
            path_t, skiprows=0, usecols=[8]) + 1j * np.loadtxt(
                path_t, skiprows=0, usecols=[11])
        Ey_t2 = np.loadtxt(
            path_t, skiprows=0, usecols=[9]) + 1j * np.loadtxt(
                path_t, skiprows=0, usecols=[12])
        Ez_t2 = np.loadtxt(
            path_t, skiprows=0, usecols=[10]) + 1j * np.loadtxt(
                path_t, skiprows=0, usecols=[13])
        Ex_t2 = Ex_t2.reshape(npt_integ, npt_integ, nb_slice, order='F')
        Ey_t2 = Ey_t2.reshape(npt_integ, npt_integ, nb_slice, order='F')
        Ez_t2 = Ez_t2.reshape(npt_integ, npt_integ, nb_slice, order='F')
        Ex_r2 = np.loadtxt(
            path_r, skiprows=0, usecols=[8]) + 1j * np.loadtxt(
                path_r, skiprows=0, usecols=[11])
        Ey_r2 = np.loadtxt(
            path_r, skiprows=0, usecols=[9]) + 1j * np.loadtxt(
                path_r, skiprows=0, usecols=[12])
        Ez_r2 = np.loadtxt(
            path_r, skiprows=0, usecols=[10]) + 1j * np.loadtxt(
                path_r, skiprows=0, usecols=[13])
        Ex_r2 = Ex_r2.reshape(npt_integ, npt_integ, nb_slice, order='F')
        Ey_r2 = Ey_r2.reshape(npt_integ, npt_integ, nb_slice, order='F')
        Ez_r2 = Ez_r2.reshape(npt_integ, npt_integ, nb_slice, order='F')

        for k11 in range(0, nb_slice):
            Ex_t3 = Ex_t2[:, :, k11]
            Ey_t3 = Ey_t2[:, :, k11]
            Ez_t3 = Ez_t2[:, :, k11]
            Ex_r3 = Ex_r2[:, :, k11]
            Ey_r3 = Ey_r2[:, :, k11]
            Ez_r3 = Ez_r2[:, :, k11]

            Ex_t3 = np.transpose(Ex_t3.conjugate())
            Ey_t3 = np.transpose(Ey_t3.conjugate())
            Ez_t3 = np.transpose(Ez_t3.conjugate())
            Ex_r3 = np.transpose(Ex_r3.conjugate())
            Ey_r3 = np.transpose(Ey_r3.conjugate())
            Ez_r3 = np.transpose(Ez_r3.conjugate())

            #plt.figure();plt.contourf(Ex_t3.real,20) ; plt.title('getdp raw Ex_t3.real (%d slice)' % k11); plt.colorbar();
            #plt.figure();plt.contourf(Ex_t3.imag,20) ; plt.title('getdp raw Ex_t3.imag (%d slice)' % k11); plt.colorbar();
            #plt.figure();plt.contourf(Ex_r3.real,20) ; plt.title('getdp raw Ex_r3.real (%d slice)' % k11); plt.colorbar();
            #plt.figure();plt.contourf(Ex_r3.imag,20) ; plt.title('getdp raw Ex_r3.imag (%d slice)' % k11); plt.colorbar();

            ex_nm_r_inter = np.zeros((1, npt_integ), dtype=complex)[0, :]
            ex_nm_t_inter = np.zeros((1, npt_integ), dtype=complex)[0, :]
            ey_nm_r_inter = np.zeros((1, npt_integ), dtype=complex)[0, :]
            ey_nm_t_inter = np.zeros((1, npt_integ), dtype=complex)[0, :]
            ez_nm_r_inter = np.zeros((1, npt_integ), dtype=complex)[0, :]
            ez_nm_t_inter = np.zeros((1, npt_integ), dtype=complex)[0, :]
            ex_nm_r = np.zeros((Nb_ordre, Nb_ordre), dtype=complex)
            ex_nm_t = np.zeros((Nb_ordre, Nb_ordre), dtype=complex)
            ey_nm_r = np.zeros((Nb_ordre, Nb_ordre), dtype=complex)
            ey_nm_t = np.zeros((Nb_ordre, Nb_ordre), dtype=complex)
            ez_nm_r = np.zeros((Nb_ordre, Nb_ordre), dtype=complex)
            ez_nm_t = np.zeros((Nb_ordre, Nb_ordre), dtype=complex)

            for n1 in range(0, Nb_ordre):
                for m1 in range(0, Nb_ordre):
                    for j1 in range(0, npt_integ):
                        expbeta = np.exp(1j * betat[m1] * y_r)
                        # ex_nm_r_inter[j1] = 1/period_y * np.trapz((Ex_r2[:,j1,k11])*expbeta,x=y_r)
                        ex_nm_r_inter[j1] = 1 / period_y * np.trapz(
                            (Ex_r3[:, j1]) * expbeta, x=y_r)
                        # plt.plot np.trapz(y_t,(Ex_t[::-1,j1].transpose()*expbeta).conjugate()[::-1])
                    expalpha = np.exp(1j * alphat[n1] * x_t)
                    ex_nm_r[n1, m1] = 1 / period_x * np.trapz(
                        ex_nm_r_inter * expalpha, x=x_r)
            for n2 in range(0, Nb_ordre):
                for m2 in range(0, Nb_ordre):
                    for j1 in range(0, npt_integ):
                        expbeta = np.exp(1j * betat[m2] * y_t)
                        # ex_nm_t_inter[j1] = 1/period_y * np.trapz((Ex_t2[:,j1,k11])*expbeta,x=y_t)
                        ex_nm_t_inter[j1] = 1 / period_y * np.trapz(
                            (Ex_t3[:, j1]) * expbeta, x=y_t)
                    expalpha = np.exp(1j * alphat[n2] * x_t)
                    ex_nm_t[n2, m2] = 1 / period_x * np.trapz(
                        ex_nm_t_inter * expalpha, x=x_t)
            for n3 in range(0, Nb_ordre):
                for m3 in range(0, Nb_ordre):
                    for j1 in range(0, npt_integ):
                        expbeta = np.exp(1j * betat[m3] * y_r)
                        # ey_nm_r_inter[j1] = 1/period_y * np.trapz((Ey_r2[:,j1,k11])*expbeta,x=y_r)
                        ey_nm_r_inter[j1] = 1 / period_y * np.trapz(
                            (Ey_r3[:, j1]) * expbeta, x=y_r)
                    expalpha = np.exp(1j * alphat[n3] * x_t)
                    ey_nm_r[n3, m3] = 1 / period_x * np.trapz(
                        ey_nm_r_inter * expalpha, x=x_r)
            for n4 in range(0, Nb_ordre):
                for m4 in range(0, Nb_ordre):
                    for j1 in range(0, npt_integ):
                        expbeta = np.exp(1j * betat[m4] * y_t)
                        # ey_nm_t_inter[j1] = 1/period_y * np.trapz((Ey_t2[:,j1,k11])*expbeta,x=y_t)
                        ey_nm_t_inter[j1] = 1 / period_y * np.trapz(
                            (Ey_t3[:, j1]) * expbeta, x=y_t)
                    expalpha = np.exp(1j * alphat[n4] * x_t)
                    ey_nm_t[n4, m4] = 1 / period_x * np.trapz(
                        ey_nm_t_inter * expalpha, x=x_t)
            for n6 in range(0, Nb_ordre):
                for m6 in range(0, Nb_ordre):
                    for j1 in range(0, npt_integ):
                        expbeta = np.exp(1j * betat[m6] * y_r)
                        # ez_nm_r_inter[j1] = 1/period_y * np.trapz((Ez_r2[:,j1,k11])*expbeta,x=y_r)
                        ez_nm_r_inter[j1] = 1 / period_y * np.trapz(
                            (Ez_r3[:, j1]) * expbeta, x=y_r)
                    expalpha = np.exp(1j * alphat[n6] * x_t)
                    ez_nm_r[n6, m6] = 1 / period_x * np.trapz(
                        ez_nm_r_inter * expalpha, x=x_r)
            for n7 in range(0, Nb_ordre):
                for m7 in range(0, Nb_ordre):
                    for j1 in range(0, npt_integ):
                        expbeta = np.exp(1j * betat[m7] * y_t)
                        # ez_nm_t_inter[j1] = 1/period_y * np.trapz((Ez_t2[:,j1,k11])*expbeta,x=y_t)
                        ez_nm_t_inter[j1] = 1 / period_y * np.trapz(
                            (Ez_t3[:, j1]) * expbeta, x=y_t)
                    expalpha = np.exp(1j * alphat[n7] * x_t)
                    ez_nm_t[n7, m7] = 1 / period_x * np.trapz(
                        ez_nm_t_inter * expalpha, x=x_t)
            ####################
            for n8 in range(0, Nb_ordre):
                for m8 in range(0, Nb_ordre):
                    AXsit[n8, m8, k11] = 1 / (layer_diopter[0]['gamma'] * gammatt[n8, m8]) * (
                        + gammatt[n8, m8]**2 * np.abs(ex_nm_t[n8, m8])**2
                        + gammatt[n8, m8]**2 * np.abs(ey_nm_t[n8, m8])**2
                        + gammatt[n8, m8]**2 * np.abs(ez_nm_t[n8, m8])**2)
            for n9 in range(0, Nb_ordre):
                for m9 in range(0, Nb_ordre):
                    AXsir[n9, m9, k11] = 1 / (layer_diopter[0]['gamma'] * gammatr[n9, m9]) * (
                        + gammatr[n9, m9]**2 * np.abs(ex_nm_r[n9, m9])**2
                        + gammatr[n9, m9]**2 * np.abs(ey_nm_r[n9, m9])**2
                        + gammatr[n9, m9]**2 * np.abs(ez_nm_r[n9, m9])**2)
        Q = self.postpro_absorption()
        Tnm = np.mean(AXsit, axis=2)
        Rnm = np.mean(AXsir, axis=2)
        # energy = dict([('trans', Tnm), ('refl', Rnm), ('abs1', Q),
        #                ('refl_slices', AXsir), ('trans_slices', AXsit)])
        balance = np.sum(np.sum(Tnm)) + np.sum(np.sum(Rnm)) + Q
        effs = dict([('T', Tnm), ('R', Rnm), ('Q', Q), ('B', balance)])
        return effs

    def ancillary_problem(self):
        nb_layer_diopter = 2
        layer_diopter = []
        for k1 in range(0, nb_layer_diopter):
            layer_diopter.append({})

        AR1 = np.zeros((1, 1), dtype=complex)[0, :]
        AT1 = np.zeros((1, 1), dtype=complex)[0, :]
        AR2 = np.zeros((1, 1), dtype=complex)[0, :]
        AT2 = np.zeros((1, 1), dtype=complex)[0, :]

        omega = 2. * pi * self.cel / self.lambda0
        k0 = 2. * pi / self.lambda0
        alpha0 = k0 * np.sin(self.theta_0) * np.cos(self.phi_0)
        beta0 = k0 * np.sin(self.theta_0) * np.sin(self.phi_0)
        gamma0 = k0 * np.cos(self.theta_0)
        gamma02 = np.sqrt(k0**2 - alpha0**2 - beta0**2)
        Ae = 1.
        Ah = Ae * np.sqrt(self.epsilon0 / self.mu0)

        Ex0 = Ae * (np.cos(self.psi_0) * np.cos(self.theta_0) * np.cos(
            self.phi_0) - np.sin(self.psi_0) * np.sin(self.phi_0))
        Ey0 = Ae * (np.cos(self.psi_0) * np.cos(self.theta_0) * np.sin(
            self.phi_0) + np.sin(self.psi_0) * np.cos(self.phi_0))
        Ez0 = Ae * (-np.cos(self.psi_0) * np.sin(self.theta_0))

        Hx0 = -1 / (omega * self.mu0) * (beta0 * Ez0 - gamma0 * Ey0)
        Hy0 = -1 / (omega * self.mu0) * (gamma0 * Ex0 - alpha0 * Ez0)
        Hz0 = -1 / (omega * self.mu0) * (alpha0 * Ey0 - beta0 * Ex0)
        #######################################################
        ####   SLAB CONFIG  (4 layer_diopters in this example)    ####
        #######################################################
        layer_diopter[0]['epsilon'] = self.eps_L1
        layer_diopter[1]['epsilon'] = self.eps_L6

        layer_diopter[0][
            'thickness'] = self.thick_L1 + self.thick_L2 + self.thick_L3 + self.thick_L4 + self.thick_L5
        layer_diopter[1]['thickness'] = self.thick_L6

        layer_diopter[0]['hh'] = 0
        layer_diopter[1]['hh'] = 0
        # for k in range(1,nb_layer_diopter):layer_diopter[k]['hh']=layer_diopter[k-1]['hh']-layer_diopter[k]['thickness']
        #################################################
        ####   SET Interface and transport matrices  ####
        #################################################
        for i_prop in range(0, nb_layer_diopter):
            layer_diopter[i_prop]['kp'] = k0 * np.sqrt(
                layer_diopter[i_prop]['epsilon'])
            layer_diopter[i_prop]['gamma'] = sc.sqrt(
                layer_diopter[i_prop]['kp']**2 - alpha0**2 - beta0**2)
            layer_diopter[i_prop]['mu'] = 1
            layer_diopter[i_prop]['M'] = sc.linalg.inv(np.array([[omega * layer_diopter[i_prop]['mu'] * self.mu0,                0,                        beta0],
                                                                 [0, omega * layer_diopter[i_prop]['mu'] *
                                                                  self.mu0,                      -alpha0],
                                                                 [-beta0,              alpha0,   -omega * layer_diopter[i_prop]['epsilon'] * self.epsilon0]]))

            layer_diopter[i_prop]['Pi'] = np.array([[1,                        1,                        0,                      0],
                                                    [0,                        0,
                                                        1,                      1],
                                                    [layer_diopter[i_prop]['gamma'] * layer_diopter[i_prop]['M'][0, 1], -layer_diopter[i_prop]['gamma'] * layer_diopter[i_prop]['M'][0, 1], -
                                                        layer_diopter[i_prop]['gamma'] * layer_diopter[i_prop]['M'][0, 0], layer_diopter[i_prop]['gamma'] * layer_diopter[i_prop]['M'][0, 0]],
                                                    [layer_diopter[i_prop]['gamma'] * layer_diopter[i_prop]['M'][1, 1], -layer_diopter[i_prop]['gamma'] * layer_diopter[i_prop]['M'][1, 1], -layer_diopter[i_prop]['gamma'] * layer_diopter[i_prop]['M'][1, 0], layer_diopter[i_prop]['gamma'] * layer_diopter[i_prop]['M'][1, 0]]])

            layer_diopter[i_prop]['T'] = np.array([[np.exp(1j * layer_diopter[i_prop]['gamma'] * layer_diopter[i_prop]['thickness']),       0,      0,       0],
                                                   [0,       np.exp(-1j * layer_diopter[i_prop]['gamma'] *
                                                                    layer_diopter[i_prop]['thickness']),      0,       0],
                                                   [0,                   0,      np.exp(
                                                       1j * layer_diopter[i_prop]['gamma'] * layer_diopter[i_prop]['thickness']),       0],
                                                   [0,                   0,                  0,       np.exp(-1j * layer_diopter[i_prop]['gamma'] * layer_diopter[i_prop]['thickness'])]])
        ##################
        ####   SOLVE  ####
        ##################
        M1 = np.eye(4)
        for i1 in range(0, nb_layer_diopter - 2):
            M1 = np.dot(
                sc.linalg.inv(layer_diopter[i1 + 1]['T']),
                np.dot(
                    sc.linalg.inv(layer_diopter[i1 + 1]['Pi']),
                    np.dot(layer_diopter[i1]['Pi'], M1)))

        M1 = np.dot(
            sc.linalg.inv(layer_diopter[nb_layer_diopter - 1]['Pi']),
            np.dot(layer_diopter[nb_layer_diopter - 2]['Pi'], M1))
        M2 = np.array([[1, 0, -M1[0, 1], -M1[0, 3]],
                       [0, 0, -M1[1, 1], -M1[1, 3]],
                       [0, 1, -M1[2, 1], -M1[2, 3]],
                       [0, 0, -M1[3, 1], -M1[3, 3]]])

        known = np.array([[M1[0, 0] * Ex0 + M1[0, 2] * Ey0],
                          [M1[1, 0] * Ex0 + M1[1, 2] * Ey0],
                          [M1[2, 0] * Ex0 + M1[2, 2] * Ey0],
                          [M1[3, 0] * Ex0 + M1[3, 2] * Ey0]])

        TetR = np.dot(sc.linalg.inv(M2), known)

        layer_diopter[nb_layer_diopter - 1]['Psi'] = np.array([TetR[0],
                                                               [0.],
                                                               TetR[1],
                                                               [0.]])

        layer_diopter[nb_layer_diopter - 2]['Psi'] = np.dot(
            sc.linalg.inv(layer_diopter[nb_layer_diopter - 2]['Pi']),
            np.dot((layer_diopter[nb_layer_diopter - 1]['Pi']),
                   layer_diopter[nb_layer_diopter - 1]['Psi']))
        for i2 in range(1, nb_layer_diopter - 1):
            layer_diopter[(nb_layer_diopter - 2) - i2]['Psi'] = np.dot(
                sc.linalg.inv(
                    layer_diopter[(nb_layer_diopter - 2) - i2]['Pi']),
                np.dot((layer_diopter[(nb_layer_diopter - 2) - i2 + 1]['Pi']),
                       np.dot((layer_diopter[(nb_layer_diopter - 2)
                                             - i2 + 1]['T']),
                              layer_diopter[(nb_layer_diopter - 2)
                                            - i2 + 1]['Psi'])))
        for i4 in range(0, nb_layer_diopter):
            layer_diopter[i4]['Psi'] = np.append(
                layer_diopter[i4]['Psi'], layer_diopter[i4]['gamma'] *
                (layer_diopter[i4]['M'][2, 0] * layer_diopter[i4]['Psi'][2] -
                 layer_diopter[i4]['M'][2, 1] * layer_diopter[i4]['Psi'][0]))
            layer_diopter[i4]['Psi'] = np.append(
                layer_diopter[i4]['Psi'], layer_diopter[i4]['gamma'] *
                (layer_diopter[i4]['M'][2, 1] * layer_diopter[i4]['Psi'][1] -
                 layer_diopter[i4]['M'][2, 0] * layer_diopter[i4]['Psi'][3]))

        AR1[0] = layer_diopter[0]['gamma'] / layer_diopter[0]['gamma'] *   \
            (abs(layer_diopter[0]['Psi'][1])**2
             + abs(layer_diopter[0]['Psi'][3])**2
             + abs(layer_diopter[0]['Psi'][5])**2)
        AT1[0] = layer_diopter[nb_layer_diopter - 1]['gamma'] / layer_diopter[0]['gamma'] *  \
            (abs(layer_diopter[nb_layer_diopter - 1]['Psi'][0])**2
             + abs(layer_diopter[nb_layer_diopter - 1]['Psi'][2])**2
             + abs(layer_diopter[nb_layer_diopter - 1]['Psi'][4])**2)

        AR2[0] = 1. / (layer_diopter[0]['gamma'] * layer_diopter[0]['gamma']) *  \
            ((layer_diopter[0]['gamma']**2 + alpha0**2) * abs(layer_diopter[0]['Psi'][1])**2
             + (layer_diopter[0]['gamma']**2 + beta0**2) *
                abs(layer_diopter[0]['Psi'][3])**2
             + 2 * alpha0 * beta0 * np.real(layer_diopter[0]['Psi'][1] * layer_diopter[0]['Psi'][3].conjugate()))
        AT2[0] = 1. / (layer_diopter[0]['gamma'] * layer_diopter[nb_layer_diopter - 1]['gamma']) *  \
            ((layer_diopter[nb_layer_diopter - 1]['gamma']**2 + alpha0**2) * abs(layer_diopter[nb_layer_diopter - 1]['Psi'][0])**2
             + (layer_diopter[nb_layer_diopter - 1]['gamma']**2 + beta0 **
                2) * abs(layer_diopter[nb_layer_diopter - 1]['Psi'][2])**2
             + 2 * alpha0 * beta0 * np.real(layer_diopter[nb_layer_diopter - 1]['Psi'][0] * layer_diopter[nb_layer_diopter - 1]['Psi'][2].conjugate()))

        # print('T_diopter=%3.7f\nR_diopter=%3.7f'%(AT2[0].real,AR2[0].real))
        # hack
        layer_diopter[0]['Trans'] = AT2[0]
        return layer_diopter

    def open_gmsh_gui(self, pos_list=["*.pos"]):
        p = [os.path.join(self.tmp_dir, pos) for pos in pos_list]
        femio.open_gmsh(
            self.path_mesh,
            self.path_geo,
            pos_list=p,
            verbose=self.gmsh_verbose)


if __name__ == "__main__":
    print("This is the femmodel module")
