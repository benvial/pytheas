import numpy as np
import scipy as sc
import os
import subprocess
import shutil
import pytheas.tools.femio as femio
from termcolor import colored
pi = np.pi


class TwoScaleFEM2D:
    """A class for the two scale convergenve homogenization of
    a 2D medium using a finite element model with Gmsh_ and GetDP_.

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

    #: flt: incident plane wave wavelength in free space
    lambda0 = 1.0

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
            if (self.getdp_verbose>=3 or self.gmsh_verbose is 4):
                sep = "-" * 51 + "\n"
            else:
                sep = ""
            toprint = sep + colored(s, 'green')
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
        param_dict["y_flag"] = int(self.y_flag)
        param_dict["save_solution"] = int(self.save_solution)
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
            des_ID, eps_des, self.content_mesh, posname, type=self.type_des)
        return femio.maketmp(eps_des_pos,  posname + ".pos", dirname=self.tmp_dir)

    def make_pos(self, des_ID, val, posname):
        ## create a pos file to be read by getdp
        self.print_progress("Creating pos file " + posname + ".pos")
        pos = femio.make_pos(
            des_ID, val, self.content_mesh, posname, type=self.type_des)
        return femio.maketmp(pos,  posname + ".pos", dirname=self.tmp_dir)


    def make_mesh(self):
        self.print_progress("Meshing model")
        femio.mesh_model(
            self.path_mesh, self.path_geo, verbose=self.gmsh_verbose)
        content_mesh = femio.get_content(self.path_mesh)
        return content_mesh

    def make_mesh_pos(self, els, nodes):
        self.print_progress("Retrieving mesh content")
        return femio.make_content_mesh_pos(nodes, els, self.dom_des,
                                           self.celltype)

    def compute_solution(self, **kwargs):
        self.print_progress("Computing solution")
        argstr = "-petsc_prealloc 1500 -ksp_type preonly \
                 -pc_type lu -pc_factor_mat_solver_package mumps"
        resolution = "electrostat_scalar"
        femio.solve_problem(
            resolution,
            self.path_pro,
            self.path_mesh,
            verbose=self.getdp_verbose,
            path_pos=self.path_pos,
            argstr=argstr)

    def ppstr(self, postop):
        return femio.postprostring(postop, self.path_pro, self.path_mesh,
                                   self.path_pos, self.getdp_verbose)

    def postpro_choice(self, name, filetype):
        if filetype in {"pos", "txt"}:
            subprocess.call(self.ppstr(name + "_" + filetype), shell=True)
        else:
            raise TypeError(
                "Wrong filetype specified: choose between txt and pos")

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

    def get_vol(self):
        return femio.load_table(self.tmp_dir + "/Vol.txt")

    def get_int_inveps(self):
        return femio.load_table(self.tmp_dir + "/I_inveps.txt")

    def postpro_fields(self, filetype="txt"):
        self.print_progress("Postprocessing fields")
        self.postpro_choice("postop_fields", filetype)



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
        fdens = sc.interpolate.NearestNDInterpolator(points,
                                                     pattern.T.flatten())
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
            eps_nodes[density == mat.threshold_val[i]] = ncomplex**2
            eps_pattern[pattern == mat.threshold_val[i]] = ncomplex**2
        return eps_nodes, eps_pattern



    def open_gmsh_gui(self, pos_list=["*.pos"]):
        self.print_progress("Opening gmsh GUI")
        p = [os.path.join(self.tmp_dir, pos) for pos in pos_list]
        femio.open_gmsh(self.path_mesh, self.path_geo, pos_list=p)



if __name__ == "__main__":
    print("This is the femmodel module")
