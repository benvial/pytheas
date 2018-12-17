import numpy as np
from pytheas import Scatt2D
import numpy.testing as npt
import os
from testutils import *


def model(verbose=False):
    fem = Scatt2D()
    fem.rm_tmp_dir()
    fem.parmesh_des = 10
    fem.parmesh = 10
    fem.parmesh_pml = 2 / 3 * fem.parmesh
    fem.matprop_pattern = [1.4, 2 - 0.02 * 1j, 3 - 0.01j]  # refractive index values
    if verbose:
        fem.getdp_verbose = 4
        fem.gmsh_verbose = 4
        fem.python_verbose = 1
    fem.initialize()
    fem.make_mesh()

    return fem


def test_scatt2D(verbose=False):
    fem = model(verbose=verbose)
    mat = pattern(xsym=True)
    fem = ref_mesh(fem, mat)
    fem.register_pattern(mat.pattern, mat._threshold_val)
    fem.compute_solution()

    ff = fem.postpro_fields_n2f()
    scs = fem.normalized_scs(ff, [0])
    scs_ref = np.array([0.267])
    # npt.assert_almost_equal(scs, scs_ref, decimal=3)

    fields_box = fem.postpro_fields_box()
    fem.postpro_fields(filetype="pos")
    fem.postpro_fields()
    filename = "u_tot.txt"
    if os.path.isfile(fem.tmppath(filename)):
        fem.get_field_map(filename)

    fields_point = fem.get_field_point()
    u_ref = -1.562936309192888 + 0.754745239167406j
    u_tot_ref = -0.5629363091928878 + 0.754745239167406j
    u_i_ref = 1 + 0j
    #
    # npt.assert_almost_equal(
    #     fields_point, np.array([u_ref, u_tot_ref, u_i_ref]), decimal=3
    # )

    # fem.adjoint = True
    #
    # fem.compute_solution()
    # fem.get_objective()
    # fem.get_adjoint()
    # fem.get_deq_deps()

    return fem
