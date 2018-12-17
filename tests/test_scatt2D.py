import numpy as np
from pytheas import Scatt2D
import numpy.testing as npt
import os
from testutils import *


def model(verbose=False):
    fem = Scatt2D()
    fem.rm_tmp_dir()
    # opto-geometric parameters  -------------------------------------------
    fem.parmesh_des = 7
    fem.parmesh = 4
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
    mat = pattern()
    fem = ref_mesh(fem, mat)
    fem.register_pattern(mat.pattern, mat._threshold_val)
    fem.compute_solution()

    ff = fem.postpro_fields_n2f()
    scs = fem.normalized_scs(ff, [0])
    scs_ref = np.array([0.19417851])
    npt.assert_almost_equal(scs, scs_ref, decimal=3)

    fields_box = fem.postpro_fields_box()
    fem.postpro_fields(filetype="pos")
    fem.postpro_fields()
    filename = "u_tot.txt"
    if os.path.isfile(fem.tmppath(filename)):
        fem.get_field_map(filename)

    fields_point = fem.get_field_point()
    u_ref = -1.615703603213303 + 0.7077065204258155j
    u_tot_ref = -0.6157036032133025 + 0.7077065204258155j
    u_i_ref = 1 + 0j

    npt.assert_almost_equal(
        fields_point, np.array([u_ref, u_tot_ref, u_i_ref]), decimal=3
    )

    fem.adjoint = True

    fem.compute_solution()
    fem.get_objective()
    fem.get_adjoint()
    fem.get_deq_deps()

    return fem
