import numpy as np
import scipy as sc
import os
import subprocess
import shutil
import pytheas.tools.femio as femio
from termcolor import colored

pi = np.pi


class BandsFEM2D:
    """A class to compute the band diagramm of a 2D
    bi-periodic medium using a finite element model with Gmsh_ and GetDP_.

        .. _Gmsh:
            http://gmsh.info/
        .. _GetDP:
            http://getdp.info/
    """

    dir_path = os.path.dirname(os.path.abspath(__file__))
    ID = "fem"  #: str: name of simulation

    # filenames
    geom_filename_ = "geometry.geo"  #: str: Gmsh geometry filename
    pro_filename_ = "main.pro"  #: str: GetDP pro filename
    bg_mesh_filename_ = "bg_mesh.geo"  #: str: Gmsh geo filename for background mesh

    content_mesh = ""
    tmp_dir = "./tmp"
    path_pos = None

    #
    # mesh_filename = "mesh.msh"  #: str: mesh filename
    # res_filename = "main.res"  #: str: GetDP res filename
    # param_filename = "parameters.dat"  #: str: Gmsh geometry filename

    getdp_verbose = 0  #: str: GetDP verbose (int between 0 and 4)
    gmsh_verbose = 0  #: str: Gmsh verbose (int between 0 and 4)
    python_verbose = 0  #: str: python verbose (int between 0 and 1)

    epsilon0 = 8.854187817e-12  #: flt: vacuum permittivity
    mu0 = 4. * pi * 1e-7  #: flt: vacuum permeability
    #: flt: speed of light in vacuum
    cel = 1.0 / (np.sqrt(epsilon0 * mu0))

    #: flt: incident wavelength in free space
    lambda0 = 1.0
    lambda_mesh = 1.0
    kx = 0  #: flt: wavevector, x component
    ky = 0  #: flt: wavevector, y component
    eps_des = 1

    #: str: polarisation (either "TE" or "TM")
    pola = "TE"

    #: flt: shift for eigenvalue search
    lambda0search = 1.0
    #: int: number of eigenvalues
    neig = 6

    #: int: number of x points for postprocessing field maps
    Nix = 51

    #: flt: global mesh parameter
    #: `MeshElementSize = lambda0/(parmesh*n)`, `n`: refractive index
    parmesh = 10.
    quad_mesh_flag = False
    type_des = "elements"
    y_flag = False
    save_solution = False

    # opto-geometric parameters  -------------------------------------------
    dx = 1  #: flt: period x
    dy = 1  #: flt: period y
    dom_des = 1000  #: design domain number (check .geo/.pro files)

    # postprocessing -------------------------------------------------
    @property
    def parmesh_des(self):
        return self.parmesh

    @property
    def Niy(self):
        return self.Nix * int(self.dy / self.dx)

    @property
    def geom_filename(self):
        return os.path.join(self.dir_path, "base", self.geom_filename_)

    @property
    def pro_filename(self):
        return os.path.join(self.dir_path, "base", self.pro_filename_)

    @property
    def bg_mesh_filename(self):
        return os.path.join(self.dir_path, "base", self.bg_mesh_filename_)

    @property
    def content_geo(self):
        return femio.get_content(self.geom_filename)

    @property
    def content_pro(self):
        return femio.get_content(self.pro_filename)

    @property
    def content_bg_mesh(self):
        return femio.get_content(self.bg_mesh_filename)

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
    def path_bg_mesh(self):
        return os.path.join(self.tmp_dir, "bg_mesh.geo")

    @property
    def path_pro(self):
        return os.path.join(self.tmp_dir, "main.pro")

    @property
    def path_mesh(self):
        return os.path.join(self.tmp_dir, "mesh.msh")

    @property
    def celltype(self):
        if self.quad_mesh_flag:
            s = "quad"
        elif not self.quad_mesh_flag:
            s = "triangle"
        return s

    @property
    def content_par(self):
        return femio.make_inputs(self.param_dict)

    def print_progress(self, s):
        if self.python_verbose:
            if self.getdp_verbose >= 3 or self.gmsh_verbose is 4:
                sep = "-" * 51 + "\n"
            else:
                sep = ""
            toprint = sep + colored(s, "green")
            print(toprint)

    def initialize(self):
        self.print_progress("Initialization")
        # tmp_name = tmp_dir.split("/")[2]
        try:
            os.mkdir(self.tmp_dir)
        except FileExistsError:
            pass
            # shutil.rmtree(tmp_dir)
        # create tmp parameters files files
        femio.maketmp(self.content_par, "parameters.dat", dirname=self.tmp_dir)
        # create tmp geo file
        femio.maketmp(self.content_geo, "geometry.geo", dirname=self.tmp_dir)
        # create tmp geo file for background mesh
        femio.maketmp(self.content_bg_mesh, "bg_mesh.geo", dirname=self.tmp_dir)
        # create tmp geo file
        femio.maketmp(self.content_pro, "main.pro", dirname=self.tmp_dir)

    def update_params(self):
        self.print_progress("Updating parameters")
        femio.maketmp(self.content_par, "parameters.dat", dirname=self.tmp_dir)

    def cleanup(self):
        os.remove("*.msh *.pre *.res *.dat *.txt *.pyc")

    #
    # def make_tmp_files(self):
    #     path_pro = femio.maketmp(pro_merged, suffix=".pro", dir=self.tmp_dir)

    def make_param_dict(self):
        param_dict = dict()
        param_dict["parmesh"] = self.parmesh
        param_dict["lambda0"] = self.lambda0
        param_dict["cel"] = self.cel
        param_dict["mu0"] = self.mu0
        param_dict["epsilon0"] = self.epsilon0
        param_dict["dx"] = self.dx
        param_dict["dy"] = self.dy
        param_dict["kx"] = self.kx
        param_dict["ky"] = self.ky
        param_dict["lambda0search"] = self.lambda0search
        param_dict["neig"] = self.neig
        param_dict["Nix"] = self.Nix
        param_dict["Niy"] = self.Niy
        return param_dict

    def generate_ID(self):
        return datetime.now().strftime("%Y_%m_%d_%H_%M_%s_%f")

    def append_ID(self, filename, extension):
        return filename + "_" + str(self.ID) + "." + extension

    def get_design_nodes(self):
        self.print_progress("Retrieving nodes")
        return femio.get_nodes(self.path_mesh, self.dom_des, self.celltype)

    def get_design_elements(self):
        self.print_progress("Retrieving elements")
        return femio.get_elements(self.path_mesh, self.dom_des, self.celltype)

    def make_eps_pos(self, des_ID, eps_des):
        ## create a pos file to be read by getdp
        posname = "eps_des"
        self.print_progress("Creating permittivity file " + posname + ".pos")
        eps_des_pos = femio.make_pos(
            des_ID, eps_des, self.content_mesh, posname, type=self.type_des
        )
        return femio.maketmp(eps_des_pos, posname + ".pos", dirname=self.tmp_dir)

    def make_pos(self, des_ID, val, posname):
        ## create a pos file to be read by getdp
        self.print_progress("Creating pos file " + posname + ".pos")
        pos = femio.make_pos(
            des_ID, val, self.content_mesh, posname, type=self.type_des
        )
        return femio.maketmp(pos, posname + ".pos", dirname=self.tmp_dir)

    def make_mesh(self):
        self.print_progress("Meshing model")
        femio.mesh_model(self.path_mesh, self.path_geo, verbose=self.gmsh_verbose)
        content_mesh = femio.get_content(self.path_mesh)
        return content_mesh

    def make_mesh_pos(self, els, nodes):
        self.print_progress("Retrieving mesh content")
        return femio.make_content_mesh_pos(nodes, els, self.dom_des, self.celltype)

    def compute_solution(self, **kwargs):
        self.print_progress("Computing solution")
        argstr = "-slepc -eps_type krylovschur \
                   -st_ksp_type preonly \
                   -st_pc_type lu \
                   -st_pc_factor_mat_solver_package mumps \
                   -eps_max_it 300 \
                   -eps_target 0.00001 \
                   -eps_target_real \
                   -eps_mpd 600 -eps_nev 400"
        resolution = self.pola
        femio.solve_problem(
            resolution,
            self.path_pro,
            self.path_mesh,
            verbose=self.getdp_verbose,
            path_pos=self.path_pos,
            argstr=argstr,
        )

    def ppstr(self, postop):
        return femio.postprostring(
            postop, self.path_pro, self.path_mesh, self.path_pos, self.getdp_verbose
        )

    def postpro_choice(self, name, filetype):
        if filetype in {"pos", "txt"}:
            subprocess.call(self.ppstr(name + "_" + filetype), shell=True)
        else:
            raise TypeError("Wrong filetype specified: choose between txt and pos")

    def postpro_fields(self, filetype="txt"):
        self.print_progress("Postprocessing fields")
        self.postpro_choice("postop_fields", filetype)

    def get_field_map(self, name):
        field = femio.load_table(self.tmp_dir + "/" + name)
        return field.reshape((self.Niy, self.Nix)).T

    def postpro_eigenvalues(self):
        self.print_progress("Retrieving eigenvalues")
        subprocess.call(self.ppstr("postop_eigenvalues_" + self.pola), shell=True)
        filename = self.tmp_dir + "/EV_" + self.pola + ".txt"
        re = np.loadtxt(filename, usecols=[1])
        im = np.loadtxt(filename, usecols=[5])
        return re + 1j * im

    def postpro_eigenvectors(self, filetype="txt"):
        self.print_progress("Retrieving eigenvectors")
        self.postpro_choice("postop_eigenvectors_" + self.pola, filetype)
        if filetype is "txt":
            filename = self.tmp_dir + "/modes_" + self.pola + ".txt"
            ev = femio.load_timetable(filename)
            return ev.reshape((self.Nix, self.Niy, self.neig))

    def make_fdens(self, pattern):
        self.print_progress("Making density function")
        n_x, n_y, n_z = pattern.shape
        x = np.linspace(-self.dx / 2, self.dx / 2, n_x + 1)
        y = np.linspace(-self.dy / 2, self.dy / 2, n_y + 1)
        z = 0
        dx, dy = x[1] - x[0], y[1] - y[0]
        x0, x1 = x[0], x[-1]
        y0, y1 = y[0], y[-1]
        x = np.linspace(x0 + dx / 2, x1 - dx / 2, n_x)
        y = np.linspace(y0 + dy / 2, y1 - dy / 2, n_y)
        xx, yy, zz = np.meshgrid(x, y, z)
        points = np.vstack((xx.ravel(), yy.ravel(), zz.ravel())).T
        fdens = sc.interpolate.NearestNDInterpolator(points, pattern.T.flatten())
        return fdens

    def assign_material(self, mat, matprop, density, lambda0):
        self.print_progress("Assigning materials")
        pattern = mat.mat_rand
        eps_nodes = np.zeros_like(density, dtype=complex)
        eps_pattern = np.zeros_like(pattern, dtype=complex)
        for i in range(mat.nb_threshold):
            if isinstance(matprop[i], str):
                ncomplex = ri.get_complex_index(lambda0, matprop[i])
            else:
                ncomplex = matprop[i]
            eps_nodes[density == mat.threshold_val[i]] = ncomplex ** 2
            eps_pattern[pattern == mat.threshold_val[i]] = ncomplex ** 2
        return eps_nodes, eps_pattern

    def open_gmsh_gui(self, pos_list=["*.pos"]):
        self.print_progress("Opening gmsh GUI")
        p = [os.path.join(self.tmp_dir, pos) for pos in pos_list]
        femio.open_gmsh(self.path_mesh, self.path_geo, pos_list=p)

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

        self.print_progress("Plotting field map")
        x = np.linspace(
            self.nper * self.domX_L, self.nper * self.domX_R, self.nper * self.Nix
        )
        y = np.linspace(self.domY_B, self.domY_T, self.Niy)
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
            alpha=0.33,
            extent=(
                self.nper * self.domX_L,
                self.nper * self.domX_R,
                self.h_layer1,
                self.h_layer1 + self.h_des,
            ),
        )
        ax.imshow(
            field,
            alpha=0.,
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

    def points_kspace(self, N):
        Gamma = [0., 0.]
        X = [1., 0.]
        M = [1., 1.]
        ngx = N
        Gamma_X = np.array(
            [np.linspace(Gamma[0], X[0], ngx), np.linspace(Gamma[1], X[1], ngx)]
        )

        nxm = N
        X_M = np.array([np.linspace(X[0], M[0], nxm), np.linspace(X[1], M[1], nxm)])

        X_M = np.delete(X_M, 0, axis=1)

        nmg = N
        M_Gamma = np.array(
            [np.linspace(M[0], Gamma[0], nmg), np.linspace(M[1], Gamma[1], nmg)]
        )
        M_Gamma = np.delete(M_Gamma, 0, axis=1)
        bandsx = np.append(np.append(Gamma_X[0, :], X_M[0, :]), M_Gamma[0, :])
        bandsy = np.append(np.append(Gamma_X[1, :], X_M[1, :]), M_Gamma[1, :])

        K1 = np.linspace(0, self.dx, ngx)

        K2 = np.linspace(self.dx, self.dx + self.dy, nxm)
        K2 = np.delete(K2, 0)
        K3 = np.linspace(
            self.dx + self.dy,
            self.dx + self.dy + np.sqrt(self.dx ** 2 + self.dy ** 2),
            nmg,
        )
        K3 = np.delete(K3, 0)
        Kplot = np.append(np.append(K1, K2), K3)
        K = np.array((bandsx, bandsy)).T
        return K, Kplot


if __name__ == "__main__":
    print("This is the femmodel module")
