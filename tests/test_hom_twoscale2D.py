import numpy as np
from pytheas import TwoScale2D

# import numpy.testing as npt
from testutils import *


def model(verbose=False):
    fem = TwoScale2D()
    fem.rm_tmp_dir()
    fem.parmesh = 50
    fem.lambda0 = fem.dx
    fem.eps_des = 1
    # opto-geometric parameters  -------------------------------------------
    fem.matprop_pattern = [1.4, 2 - 0.02 * 1j, 3 - 0.01j]  # refractive index values
    if verbose:
        fem.getdp_verbose = 4
        fem.gmsh_verbose = 4
        fem.python_verbose = 1
    fem.initialize()
    fem.make_mesh()

    return fem


def test_hom(verbose=False):
    fem = model(verbose=verbose)
    mat = pattern(xsym=False, ysym=False)
    # fem = ref_mesh(fem, mat, periodic_x=True, periodic_y=True)
    fem.register_pattern(mat.pattern, mat._threshold_val)
    fem.compute_solution()
    fem.compute_epsilon_eff()
    # eps_eff_ref = np.array(
    #     [
    #         [4.09598751 - 0.06423691j, 0.0353099 - 0.001815j],
    #         [0.0353099 - 0.001815j, 3.82597454 - 0.05819867j],
    #     ]
    # )

    # npt.assert_almost_equal(eps_eff, eps_eff_ref, decimal=2)
    return fem
