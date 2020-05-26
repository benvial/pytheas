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

import shutil
import os
import subprocess
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from ..tools import femio


def get_file_path(f):
    return os.path.dirname(os.path.abspath(f))


class BaseFEM:
    """Base class for Finite Element models
    """

    epsilon0 = 8.854187817e-12  #: flt: vacuum permittivity
    mu0 = 4.0 * np.pi * 1e-7  #: flt: vacuum permeability
    cel = 1.0 / (np.sqrt(epsilon0 * mu0))  #: flt: speed of light in vacuum

    def __init__(self):
        self.geom_filename_ = "geometry.geo"  #: str: Gmsh geometry filename
        self.pro_filename_ = "main.pro"  #: str: GetDP pro filename
        self.param_filename_ = "parameters.dat"  #: str: GetDP pro filename
        #: str: Gmsh geo filename for background mesh
        self.bg_mesh_filename_ = "bg_mesh.geo"
        self.bg_mesh = True
        # : bool: wether or not to use an inclusion geometry instead of a material distribution
        self.inclusion_flag = False
        self.pola = None
        self.adjoint = False

        self.inclusion_filename_ = "inclusion.geo"
        self.content_mesh = ""
        self.tmp_dir = "./tmp"
        self.path_pos = ""
        self.getdp_verbose = 0  #: str: GetDP verbose (int between 0 and 4)
        self.gmsh_verbose = 0  #: str: Gmsh verbose (int between 0 and 4)
        self.python_verbose = 0  #: str: python verbose (int between 0 and 1)
        #: flt: global mesh parameter
        #: `MeshElementSize = lambda0/(parmesh*n)` `n`: refractive index
        self.parmesh = 10.0
        #: flt: design subdomain mesh parameter
        self.parmesh_des = 10.0
        #: flt: PMLs mesh parameter
        self.parmesh_pml = 7.0
        self.parmesh_incl = 10.0
        self.dim = 2  #: dimension of the problem
        self.quad_mesh_flag = False
        self.extrude_mesh_flag = False
        self.type_des = "elements"
        #: int: number of x points for postprocessing field maps
        self.Nix = 100
        self.Niy = 100
        self.matprop_pattern = 0
        self.pattern = False
        self.cplx_list = ["eps_"]
        self.dom_des = 0
        self.param_dict = dict()
        self.dir_path = get_file_path(__file__)
        self.analysis = "direct"
        self.ignore_periodicity = False
        self._debug = False
        self.npt_integ = 200

    @property
    def debug(self):
        return self._debug

    @debug.setter
    def debug(self, val):
        self._debug = val
        if self._debug:
            self.getdp_verbose = 4
            self.gmsh_verbose = 4
            self.python_verbose = 1
        return self._debug

    @property
    def geom_filename(self):
        return os.path.join(self.dir_path, "base", self.geom_filename_)

    @property
    def inclusion_filename(self):
        return self.tmppath(self.inclusion_filename_)

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
        if self.bg_mesh:
            return femio.get_content(self.bg_mesh_filename)
        else:
            return ""

    #
    # @property
    # def param_dict(self):
    #     return self._make_param_dict()

    # @property
    # def dir_path(self):
    #     return os.path.dirname(os.path.abspath(__file__))

    @property
    def content_par(self):
        return femio.make_inputs(self.param_dict)

    @property
    def path_geo(self):
        return self.tmppath(self.geom_filename_)

    @property
    def path_bg_mesh(self):
        return self.tmppath(self.bg_mesh_filename_)

    @property
    def path_pro(self):
        return self.tmppath(self.pro_filename_)

    @property
    def path_mesh(self):
        return self.tmppath("mesh.msh")

    @property
    def celltype(self):
        if self.quad_mesh_flag:
            s = "quad"
        elif not self.quad_mesh_flag:
            s = "triangle"
        return s

    def tmppath(self, f):
        return os.path.join(self.tmp_dir, f)

    def _print_progress(self, s):
        if self.python_verbose:
            if self.getdp_verbose >= 3 or self.gmsh_verbose is 4:
                sep = "-" * 51 + "\n"
            else:
                sep = ""
            print(sep + s)

    def initialize(self):
        """Initialize the problem parameters.
        """
        self._print_progress("Initialization")
        # tmp_name = tmp_dir.split("/")[2]
        self.mk_tmp_dir()

        # create tmp parameters files files
        self.param_dict = self._make_param_dict()
        femio.maketmp(self.content_par, self.param_filename_, dirname=self.tmp_dir)
        # create tmp geo file
        femio.maketmp(self.content_geo, self.geom_filename_, dirname=self.tmp_dir)
        # create tmp geo file for background mesh
        if self.bg_mesh:
            femio.maketmp(
                self.content_bg_mesh, self.bg_mesh_filename_, dirname=self.tmp_dir
            )
        # create tmp pro file
        femio.maketmp(self.content_pro, self.pro_filename_, dirname=self.tmp_dir)
        # if self.inclusion_flag:
        #     # create tmp geo inclusion file
        #     femio.maketmp(self.content_incl, "inclusion.geo", dirname=self.tmp_dir)

    def update_params(self):
        """
        Update the dictionary of parameters and the corresponding file
        """
        self._print_progress("Updating parameters")
        self.param_dict = self._make_param_dict()
        femio.maketmp(self.content_par, self.param_filename_, dirname=self.tmp_dir)

    def cleanup(self):
        """Remove gmsh/getdp/python generated files from the temporary folder"""
        trash = ["*.msh", "*.pre", "*.res", "*.dat", "*.txt", "*.pyc", "*.pos"]
        for item in trash:
            try:
                os.remove(self.tmppath(item))
            except OSError:
                pass
        return

    def mk_tmp_dir(self):
        """Create a temporary directory"""
        try:
            os.mkdir(self.tmp_dir)
            if self.python_verbose:
                print("Creating temporary directory {}".format(self.tmp_dir))
        except FileExistsError as er:
            if self.python_verbose:
                print(er)
                print("Writing inside...")
            else:
                pass
        return

    def rm_tmp_dir(self):
        """Remove the temporary directory"""
        try:
            shutil.rmtree(self.tmp_dir)
            if self.python_verbose:
                print("Removed temporary directory {}".format(self.tmp_dir))
        except FileNotFoundError as er:
            if self.python_verbose:
                print(er)
                print("Skipping...")
            else:
                pass
        return

    def _make_param_dict(self):
        """Build dictionary of parameters. This will be later written to a parameter.dat
        file that is meant to be read by both gmsh and getdp"""
        param_dict = dict()
        attr_list = [i for i in dir(self) if i[:1] != "_"]
        attr_list = [i for i in attr_list if not callable(getattr(self, i))]
        for key, val in self.__dict__.items():
            for cpl in self.cplx_list:
                if key.startswith(cpl):
                    if isinstance(val, (float, np.float64, int)):
                        self.__dict__[key] = complex(val)
        for key in attr_list:
            val = getattr(self, key)
            if isinstance(val, complex):
                param_dict[key + "_re"] = val.real
                param_dict[key + "_im"] = val.imag
            if isinstance(val, bool):
                # special handling
                param_dict[key] = int(val)
                pass
            elif isinstance(val, (float, int)):
                param_dict[key] = val
            else:
                pass

        param_dict["TE_flag"] = int(self.pola == "TE")
        param_dict["inclusion_flag"] = int(self.inclusion_flag)
        param_dict["adjoint_flag"] = int(self.adjoint)
        param_dict["quad_mesh_flag"] = int(self.quad_mesh_flag)
        param_dict["extrude_mesh_flag"] = int(self.extrude_mesh_flag)
        param_dict["nodes_flag"] = int(self.type_des == "nodes")
        return param_dict

    def make_inclusion(self, points, lcar="lc_incl", **kwargs):
        """Make a diffractive element geometry from points.

        Parameters
        ----------
        points : array of size (Npoints, 2)
            The points defining the simply connected 2D geometry of the object.
        lcar : str (default "lc_incl")
            Caracteristic length for the mesh.
        **kwargs : dict
            Extra arguments.

        """
        femio.points2geo(points, lcar, output_path=self.inclusion_filename, **kwargs)

    def get_design_nodes(self):
        self._print_progress("Retrieving nodes")
        return femio.get_nodes(self.path_mesh, self.dom_des, self.celltype)

    def get_design_elements(self):
        self._print_progress("Retrieving elements")
        return femio.get_elements(self.path_mesh, self.dom_des, self.celltype)

    def make_eps_pos(self, des_ID, _eps_des, posname="eps_des"):
        # create a pos file to be read by getdp
        self._print_progress("Creating permittivity file " + posname + ".pos")
        eps_des_pos = femio.make_pos(
            des_ID, _eps_des, self.content_mesh, posname, celltype=self.type_des
        )
        return femio.maketmp(eps_des_pos, posname + ".pos", dirname=self.tmp_dir)

    def make_pos(self, des_ID, val, posname):
        # create a pos file to be read by getdp
        self._print_progress("Creating pos file " + posname + ".pos")
        pos = femio.make_pos(
            des_ID, val, self.content_mesh, posname, celltype=self.type_des
        )
        return femio.maketmp(pos, posname + ".pos", dirname=self.tmp_dir)

    def make_mesh(self, other_option=None):
        """Mesh the geometry using gmsh.

        Parameters
        ----------
        other_option : str
            Extra flag to pass to gmsh.

        Returns
        -------
        str
            The content of the .msh file.

        """
        other_option = other_option or ""

        if self.dim == 3:
            dim = [1, 2, 3]
        else:
            dim = [1, 2]
        self._print_progress("Meshing model")
        if self.ignore_periodicity:
            print("Ignoring periodicity")
            igper = "-ignore_periodicity"
        else:
            igper = ""
        femio.mesh_model(
            self.path_mesh,
            self.path_geo,
            dim=dim,
            verbose=self.gmsh_verbose,
            other_option=other_option + igper,
        )
        self.content_mesh = femio.get_content(self.path_mesh)
        self.get_mesh_info()
        return self.content_mesh

    def make_mesh_pos(self, els, nodes):
        self._print_progress("Retrieving mesh content")
        return femio.make_content_mesh_pos(nodes, els, self.dom_des, self.celltype)

    def compute_solution(self, res_list=None):
        """Compute the solution of the FEM problem using getdp"""
        res_list = res_list or ["helmoltz_scalar", "helmoltz_scalar_modal"]
        if self.pattern:
            self.update_epsilon_value()
        self.update_params()
        self._print_progress("Computing solution: " + self.analysis + " problem")
        if self.analysis == "direct":
            argstr = "-petsc_prealloc 200 -ksp_type preonly \
                     -pc_type lu -pc_factor_mat_solver_type mumps"

            resolution = res_list[0]
        elif self.analysis == "modal":
            argstr = "-slepc -eps_type krylovschur \
                       -st_ksp_type preonly \
                       -st_pc_type lu \
                       -st_pc_factor_mat_solver_type mumps \
                       -eps_max_it 300 \
                       -eps_target 0.00001 \
                       -eps_target_real \
                       -eps_mpd 600 -eps_nev 400"

            resolution = res_list[1]
        else:
            raise TypeError("Wrong analysis specified: choose between direct and modal")

        argstr += " -cpu"
        femio.solve_problem(
            resolution,
            self.path_pro,
            self.path_mesh,
            verbose=self.getdp_verbose,
            path_pos=self.path_pos,
            argstr=argstr,
        )

    def _ppcmd(self, postop):
        """Create a postprocessing command

        Parameters
        ----------
        postop : str
            Name of the post operation as defined in the .pro file.
        """
        return femio.postpro_commands(
            postop, self.path_pro, self.path_mesh, self.path_pos, self.getdp_verbose
        )

    def _postpro_choice(self, name, filetype):
        """Run a postprocessing command with either 'pos' or 'txt' file output.

        Parameters
        ----------
        name : str
            Name of the post operation as defined in the .pro file.
        filetype : str
            File type to use ('pos' or 'txt')
        """

        if filetype in {"pos", "txt"}:
            self.postprocess(name + "_" + filetype)
        else:
            raise TypeError("Wrong filetype specified: choose between txt and pos")

    def _get_qty(self, filename):
        """Retrieve a scalar quantity.

        Parameters
        ----------
        filename : str
            Name of the txt file to load.

        Returns
        -------
        qty : array
            The quantity to be loaded.
        """
        file_path = self.tmppath(filename)
        if self.type_des is "nodes":
            return femio.load_node_table(file_path)[1]
        else:
            return femio.load_table(file_path)

    def get_qty_vect(self, filename):
        file_path = self.tmppath(filename)
        if self.type_des is "nodes":
            return femio.load_node_table_vect(file_path)[1]
        else:
            return femio.load_table_vect(file_path)

    def make_fdens(self, pattern):
        self._print_progress("Making density function")
        n_x, n_y, n_z = pattern.shape
        if len(self.corners_des) == 6:
            x0, x1, y0, y1, z0, z1 = self.corners_des
        else:
            x0, x1, y0, y1 = self.corners_des
        x = np.linspace(x0, x1, n_x + 1)
        y = np.linspace(y0, y1, n_y + 1)
        dx, dy = x[1] - x[0], y[1] - y[0]
        if len(self.corners_des) == 6:
            z = np.linspace(z0, z1, n_z + 1)
            dz = z[1] - z[0]
        else:
            z0, z1 = 0, 0
            dz = 0
        x = np.linspace(x0 + dx / 2, x1 - dx / 2, n_x)
        y = np.linspace(y0 + dy / 2, y1 - dy / 2, n_y)
        z = np.linspace(z0 + dz / 2, z1 - dz / 2, n_z)
        xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
        points = np.vstack((xx.ravel(), yy.ravel(), zz.ravel())).T
        fdens = NearestNDInterpolator(points, pattern.flatten())
        return fdens

    def assign_material(self, mat, matprop, density, lambda0):
        self._print_progress("Assigning materials")
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

    def get_mesh_info(self):
        # get nodes and elements and their IDs in the design domain
        nodes = self.get_design_nodes()
        els = self.get_design_elements()
        nodes_ID, nodes_coords = nodes
        els_ID, els_coords, _, _ = els
        if self.type_des is "elements":
            des = els_ID, els_coords
        elif self.type_des is "nodes":
            des = nodes_ID, nodes_coords
        self.nodes, self.els, self.des = nodes, els, des
        return nodes, els, des

    def register_pattern(self, pattern, threshold_val):
        self.pattern_ = pattern
        self.threshold_val = threshold_val
        # self.content_mesh = self.make_mesh_pos(self.els, self.nodes)
        # define a density function from a pattern
        self.fdens = self.make_fdens(pattern)
        # interpolate
        self.density = self.fdens(self.des[1])
        self.pattern = True

    def update_epsilon_value(self):
        self._print_progress("Assigning materials")
        # assign the permittivity
        self._eps_des, self.eps_pattern = assign_epsilon(
            self.pattern_, self.matprop_pattern, self.threshold_val, self.density
        )
        # create a pos file to be read by getdp
        self.path_pos = self.make_eps_pos(self.des[0], self._eps_des)

    def open_gmsh_gui(self, pos_list=None):
        """Open gmsh GUI to visualize geometry and postprocessing results.

        Parameters
        ----------
        pos_list : list
            A list of .pos files giving the views to load. By default it will
            render all the generated views.

        """
        pos_list = pos_list or ["*.pos"]
        self._print_progress("Opening gmsh GUI")
        p = [self.tmppath(pos) for pos in pos_list]
        femio.open_gmsh(self.path_mesh, self.path_geo, pos_list=p)

    def postpro_eigenvalues(
        self, postop="postop_eigenvalues", eig_file="EigenValues.txt"
    ):
        self._print_progress("Retrieving eigenvalues")
        self.postprocess(postop)
        filename = self.tmppath(eig_file)
        return femio.load_ev_timetable(filename)

    def postpro_eigenvectors(
        self, filetype="txt", postop="postop_eigenvectors", eig_file="EigenVectors.txt"
    ):
        self._print_progress("Retrieving eigenvectors")
        self._postpro_choice(postop, filetype)
        if filetype is "txt":
            filename = self.tmppath(eig_file)
            mode = femio.load_timetable(filename)
            u1 = np.zeros((self.Nix, self.Niy, self.neig), dtype=complex)
            u = mode.reshape((self.Niy, self.Nix, self.neig))
            for imode in range(self.neig):
                u1[:, :, imode] = np.flipud(u[:, :, imode]).T
            return u1
        else:
            return

    def get_spectral_elements(self):
        eigval = self.postpro_eigenvalues()
        eigvect = self.postpro_eigenvectors()
        isort = np.argsort(eigval)
        eigval = eigval[isort]
        eigvect = eigvect[:, :, (isort)]
        return eigval, eigvect

    def postpro_norm_eigenvectors(
        self, postop="postop_norm_eigenvectors", eig_file="NormsEigenVectors.txt"
    ):
        self._print_progress("Retrieving eigenvector norms")
        self.postprocess(postop)
        filename = self.tmppath(eig_file)
        return np.sqrt(femio.load_timetable(filename))

    def postprocess(self, postop):
        """Run getdp postoperation.

        Parameters
        ----------
        postop : str
            Name of the postoperation to run.
        """
        subprocess.call(self._ppcmd(postop))

    def postpro_fields(self, filetype="txt", postop="postop_fields"):
        """ Compute the field maps and output to a file.

            Parameters
            ----------
            filetype : str, default "txt"
                Type of output files. Either "txt" (to be read by the method
                get_field_map in python) or "pos" to be read by gmsh/getdp.
            postop : str, default "postop_fields"
                Name of the postoperation

        """
        self._print_progress("Postprocessing fields")
        self._postpro_choice(postop, filetype)

    def postpro_fields_pos(self, postop="postop_fields"):
        return self.postpro_fields(filetype="pos", postop=postop)

    def get_objective(self, postop="postop_int_objective", filename="objective.txt"):
        self._print_progress("Retrieving objective")
        if not self.adjoint:
            self.postprocess(postop)
        return femio.load_table(self.tmppath(filename)).real

    def get_adjoint(self, name="adjoint.txt"):
        self._print_progress("Retrieving adjoint")
        if self.dim is 2:
            return self._get_qty(name)
        else:
            return self.get_qty_vect(name)

    def get_deq_deps(self, name="dEq_deps.txt"):
        self._print_progress("Retrieving dEq_deps")
        if self.dim is 2:
            return self._get_qty(name)
        else:
            return self.get_qty_vect(name)


def assign_epsilon(pattern, matprop, threshold_val, density):
    _eps_des = np.zeros_like(density, dtype=complex)
    eps_pattern = np.zeros_like(pattern, dtype=complex)
    for ncomplex, thres in zip(matprop, threshold_val):
        _eps_des[density == thres] = ncomplex ** 2
        eps_pattern[pattern == thres] = ncomplex ** 2
    return _eps_des, eps_pattern
