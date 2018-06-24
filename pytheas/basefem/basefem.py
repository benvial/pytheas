# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT

"""
Base class for Finite Element models
====================================

Define, solve and postprocess a FEM model using Gmsh_ and GetDP_.

 .. _Gmsh:
     http://gmsh.info/
 .. _GetDP:
     http://getdp.info/
"""


import os
import subprocess
import numpy as np
import scipy as sc
import tempfile
from termcolor import colored
from ..tools import femio

# import sys
# sys.path.insert(0, '/home/bench/Codes/python/packages/custom_inherit')
# from custom_inherit import DocInheritMeta
#
# class BaseFEM(metaclass=DocInheritMeta(style="numpy")):
class BaseFEM():
    """Base class for Finite Element models

    Parameters
    ----------
    ID : str, default 'sim'
        Name of the simulation
    """

    epsilon0 = 8.854187817e-12  #: flt: vacuum permittivity
    mu0 = 4. * np.pi * 1e-7  #: flt: vacuum permeability
    cel = 1.0 / (np.sqrt(epsilon0 * mu0))  #: flt: speed of light in vacuum

    dir_path = os.path.dirname(os.path.abspath(__file__))

    def __init__(self, ID="sim"):
        self.ID = ID
        self.geom_filename_ = "geometry.geo"  #: str: Gmsh geometry filename
        self.pro_filename_ = "main.pro"  #: str: GetDP pro filename
        #: str: Gmsh geo filename for background mesh
        self.bg_mesh_filename_ = "bg_mesh.geo"
        # : bool: wether or not to use an inclusion geometry instead of a material distribution
        self.inclusion_flag = False

        self.inclusion_filename_ = "inclusion.geo"
        self.content_mesh = ""
        self.tmp_dir = "./tmp"
        self.path_pos = None
        self.getdp_verbose = 0  #: str: GetDP verbose (int between 0 and 4)
        self.gmsh_verbose = 0  #: str: Gmsh verbose (int between 0 and 4)
        self.python_verbose = 0  #: str: python verbose (int between 0 and 1)
        #: flt: global mesh parameter
        #: `MeshElementSize = lambda0/(parmesh*n)` `n`: refractive index
        self.parmesh = 10.
        #: flt: design subdomain mesh parameter
        self.parmesh_des = 10.
        #: flt: PMLs mesh parameter
        self.parmesh_pml = 7.
        self.parmesh_incl = 10.
        self.quad_mesh_flag = False
        self.type_des = "elements"
        #: int: number of x points for postprocessing field maps
        self.Nix = 100
        self.Niy = 100
        self.matprop_pattern = 0
        self.pattern = False
        self.param_dict = dict()

    @property
    def geom_filename(self):
        return os.path.join(self.dir_path, "base", self.geom_filename_)

    @property
    def inclusion_filename(self):
        return os.path.join(self.tmp_dir, self.inclusion_filename_)

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

    #
    # @property
    # def param_dict(self):
    #     return self.make_param_dict()

    @property
    def content_par(self):
        return femio.make_inputs(self.param_dict)

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

    def print_progress(self, s):
        if self.python_verbose:
            if (self.getdp_verbose >= 3 or self.gmsh_verbose is 4):
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
        self.param_dict = self.make_param_dict()
        femio.maketmp(self.content_par, "parameters.dat", dirname=self.tmp_dir)
        # create tmp geo file
        femio.maketmp(self.content_geo, "geometry.geo", dirname=self.tmp_dir)
        # create tmp geo file for background mesh
        femio.maketmp(self.content_bg_mesh,
                      "bg_mesh.geo", dirname=self.tmp_dir)
        # create tmp pro file
        femio.maketmp(self.content_pro, "main.pro", dirname=self.tmp_dir)
        # if self.inclusion_flag:
        #     # create tmp geo inclusion file
        #     femio.maketmp(self.content_incl, "inclusion.geo", dirname=self.tmp_dir)

    def update_params(self):
        self.print_progress("Updating parameters")
        self.param_dict = self.make_param_dict()
        femio.maketmp(self.content_par, "parameters.dat", dirname=self.tmp_dir)

    def cleanup(self):
        """Clean gmsh/getdp generated files

        """
        trash = ["*.msh", "*.pre", "*.res", "*.dat", "*.txt", "*.pyc", "*.pos"]
        for item in trash:
            try:
                os.remove(os.path.join(self.tmp_dir, item))
            except OSError:
                pass

    def make_param_dict(self):
        param_dict = dict()
        attr_list = [i for i in dir(self) if i[:1] != '_']
        attr_list = [i for i in attr_list if not callable(getattr(self, i))]
        for key, val in self.__dict__.items():
            if key.startswith("eps_"):
                if type(val) is float or type(val) is np.float64 or type(val) is int:
                    self.__dict__[key] = complex(val)
        for key in attr_list:
            val = getattr(self, key)
            if isinstance(val, complex):
                param_dict[key + "_re"] = val.real
                param_dict[key + "_im"] = val.imag
            if isinstance(val, bool):
                # special handling
                pass
            elif isinstance(val, (float, int)):
                param_dict[key] = val
            else:
                pass

        param_dict["TE_flag"] = int(self.pola == "TE")
        param_dict["inclusion_flag"] = int(self.inclusion_flag)
        param_dict["adjoint_flag"] = int(self.adjoint)
        param_dict["quad_mesh_flag"] = int(self.quad_mesh_flag)
        param_dict["nodes_flag"] = int(self.type_des == "nodes")
        return param_dict

    def generate_ID(self):
        return datetime.now().strftime("%Y_%m_%d_%H_%M_%s_%f")

    def append_ID(self, filename, extension):
        return filename + "_" + str(self.ID) + "." + extension

    def make_inclusion(self, points, **kwargs):
        femio.points2geo(points, "lc_incl",
                         output_path=self.inclusion_filename, **kwargs)

    def get_design_nodes(self):
        self.print_progress("Retrieving nodes")
        return femio.get_nodes(self.path_mesh, self.dom_des, self.celltype)

    def get_design_elements(self):
        self.print_progress("Retrieving elements")
        return femio.get_elements(self.path_mesh, self.dom_des, self.celltype)

    def make_eps_pos(self, des_ID, _eps_des):
        # create a pos file to be read by getdp
        posname = "eps_des"
        self.print_progress("Creating permittivity file " + posname + ".pos")
        eps_des_pos = femio.make_pos(
            des_ID, _eps_des, self.content_mesh, posname, type=self.type_des)
        return femio.maketmp(eps_des_pos,  posname + ".pos", dirname=self.tmp_dir)

    def make_pos(self, des_ID, val, posname):
        # create a pos file to be read by getdp
        self.print_progress("Creating pos file " + posname + ".pos")
        pos = femio.make_pos(
            des_ID, val, self.content_mesh, posname, type=self.type_des)
        return femio.maketmp(pos,  posname + ".pos", dirname=self.tmp_dir)

    def make_mesh(self, other_option=" -cpu"):
        self.print_progress("Meshing model")
        femio.mesh_model(
            self.path_mesh, self.path_geo, verbose=self.gmsh_verbose,
            other_option=other_option)
        self.content_mesh = femio.get_content(self.path_mesh)
        self.nodes, self.els, self.des = self.get_mesh_info()
        return self.content_mesh


    def make_mesh_pos(self, els, nodes):
        self.print_progress("Retrieving mesh content")
        return femio.make_content_mesh_pos(nodes, els, self.dom_des,
                                           self.celltype)

    def compute_solution(self, **kwargs):
        if self.pattern:
            self.update_epsilon_value()
        self.update_params()
        self.print_progress("Computing solution: " +
                            self.analysis + " problem")
        if self.analysis == "diffraction":
            argstr = "-petsc_prealloc 1500 -ksp_type preonly \
                     -pc_type lu -pc_factor_mat_solver_package mumps"

            resolution = "helmoltz_scalar"
        elif self.analysis == "modal":
            argstr = "-slepc -eps_type krylovschur \
                       -st_ksp_type preonly \
                       -st_pc_type lu \
                       -st_pc_factor_mat_solver_package mumps \
                       -eps_max_it 300 \
                       -eps_target 0.00001 \
                       -eps_target_real \
                       -eps_mpd 600 -eps_nev 400"

            resolution = "helmoltz_scalar_modal"
        elif self.analysis == "electrostatic":
            argstr = "-petsc_prealloc 1500 -ksp_type preonly \
                     -pc_type lu -pc_factor_mat_solver_package mumps"

            resolution = "electrostat"
        else:
            raise TypeError(
                "Wrong analysis specified: choose between diffraction, modal and electrostatic"
            )

        argstr += " -cpu"
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

    def get_qty(self, filename):
        file_path = os.path.join(self.tmp_dir, filename)
        if self.type_des is "nodes":
            a = femio.load_node_table(file_path)[1]
        else:
            a = femio.load_table(file_path)
        return a

    def make_fdens(self, pattern):
        self.print_progress("Making density function")
        n_x, n_y, n_z = pattern.shape
        x0, x1, y0, y1 = self.corners_des
        x = np.linspace(x0, x1, n_x + 1)
        y = np.linspace(y0, y1, n_y + 1)
        z = 0
        dx, dy = x[1] - x[0], y[1] - y[0]
        x = np.linspace(x0 + dx / 2, x1 - dx / 2, n_x)
        y = np.linspace(y0 + dy / 2, y1 - dy / 2, n_y)
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        points = np.vstack((xx.ravel(), yy.ravel(), zz.ravel())).T
        fdens = sc.interpolate.NearestNDInterpolator(points,
                                                     pattern.flatten())
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


    def get_mesh_info(self):
        # get nodes and elements and their IDs in the design domain
        nodes = self.get_design_nodes()
        els = self.get_design_elements()
        nodes_ID, nodes_coords = nodes
        els_ID, els_coords, els_nodes_ID, geom_ID_dom = els
        xnodes, ynodes, znodes = nodes_coords.T
        xels, yels, zels = els_coords.T
        if self.type_des is "elements":
            des_ID, des_coords = els_ID, els_coords
            des = els_ID, els_coords
        elif self.type_des is "nodes":
            des_ID, des_coords = nodes_ID, nodes_coords
            des = nodes_ID, nodes_coords
        return nodes, els, des



    def register_pattern(self, pattern, threshold_val):
        self.pattern = pattern
        self.threshold_val = threshold_val
        self.content_mesh = self.make_mesh_pos(self.els, self.nodes)
        # define a density function from a pattern
        self.fdens = self.make_fdens(pattern)
        # interpolate
        self.density = self.fdens(self.des[1])
        self.pattern = True

    def update_epsilon_value(self):
        self.print_progress("Assigning materials")
        # assign the permittivity
        self._eps_des, self.eps_pattern = assign_epsilon(self.pattern, self.matprop_pattern,
                                                   self.threshold_val, self.density)
        # create a pos file to be read by getdp
        self.path_pos = self.make_eps_pos(self.des[0], self._eps_des)


    def open_gmsh_gui(self, pos_list=["*.pos"]):
        self.print_progress("Opening gmsh GUI")
        p = [os.path.join(self.tmp_dir, pos) for pos in pos_list]
        femio.open_gmsh(self.path_mesh, self.path_geo, pos_list=p)



def assign_epsilon(pattern, matprop, threshold_val, density):
    _eps_des = np.zeros_like(density, dtype=complex)
    eps_pattern = np.zeros_like(pattern, dtype=complex)
    for ncomplex, thres in zip(matprop, threshold_val):
        _eps_des[density == thres] = ncomplex**2
        eps_pattern[pattern == thres] = ncomplex**2
    return _eps_des, eps_pattern