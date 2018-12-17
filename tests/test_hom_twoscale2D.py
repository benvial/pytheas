import numpy as np
from pytheas import TwoScale2D
import numpy.testing as npt
from testutils import *


def model(verbose=False):
    fem = TwoScale2D()
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
    mat = pattern()
    fem.register_pattern(mat.pattern, mat._threshold_val)
    fem.compute_solution()
    eps_eff = fem.compute_epsilon_eff()
    eps_eff_ref = np.array(
        [
            [5.91250942e+00 - 0.07246627j, 5.97109657e-03 + 0.00041939j],
            [5.97109657e-03 + 0.00041939j, 6.10241637e+00 - 0.06557491j],
        ]
    )
    npt.assert_almost_equal(eps_eff, eps_eff_ref, decimal=3)
    return fem
