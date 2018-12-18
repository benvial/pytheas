import numpy as np
from pytheas import Scatt3D

# import numpy.testing as npt
from testutils import *


def model(verbose=False):
    fem = Scatt3D()
    fem.rm_tmp_dir()
    fem.eps_des = 2 - 0 * 1j
    fem.hx_des = 3
    fem.hy_des = 3
    fem.hz_des = 0.5
    fem.hx_box = fem.hx_des * 1.1
    fem.hy_box = fem.hy_des * 1.1
    fem.hz_box = fem.hz_des * 3
    fem.R_sph = 0.05
    fem.z_sph = -fem.hz_des / 2 * 1.2 - fem.R_sph
    fem.parmesh = 2
    fem.parmesh_des = 2
    fem.parmesh_pml = fem.parmesh * 2 / 3
    fem.parmesh_host = 3
    fem.recomb_des = False
    fem.matprop_pattern = [1.4, 2 - 0.02 * 1j, 3 - 0.01j]  # refractive index values
    if verbose:
        fem.getdp_verbose = 4
        fem.gmsh_verbose = 4
        fem.python_verbose = 1
    fem.initialize()
    fem.make_mesh()

    return fem


def test_scatt3D(verbose=False):
    fem = model(verbose=verbose)
    mat = pattern(xsym=True, ysym=True)
    fem.register_pattern(mat.pattern, mat._threshold_val)
    fem.compute_solution()
    return fem
