import numpy as np
from pytheas import BandDiag2D
import numpy.testing as npt
from testutils import *
import matplotlib.pyplot as plt


def model(verbose=False):
    fem = BandDiag2D()
    # opto-geometric parameters  -------------------------------------------
    fem.dx = 1
    fem.dy = 1
    fem.lambda0search = 10
    fem.neig = 11
    fem.Nix = 101
    fem.parmesh = 10
    fem.type_des = "elements"
    fem.matprop_pattern = [1.4, 2 - 0.02 * 1j, 3 - 0.01j]  # refractive index values
    if verbose:
        fem.getdp_verbose = 4
        fem.gmsh_verbose = 4
        fem.python_verbose = 1
    fem.initialize()
    fem.make_mesh()
    return fem


def test_eigpb(verbose=False):
    fem = model(verbose=verbose)
    mat = pattern()
    fem.register_pattern(mat.pattern, mat._threshold_val)
    fem.kx, fem.ky = 0, 0
    fem.pola = "TE"
    # fem.initialize()
    fem.compute_solution()
    _ = fem.postpro_eigenvalues()
    evTE = np.sort(_) * fem.dx / (2 * np.pi)
    fem.pola = "TM"
    fem.update_params()
    fem.compute_solution()
    _ = fem.postpro_eigenvalues()
    evTM = np.sort(_) * fem.dx / (2 * np.pi)

    evTE_ref = np.array(
        [
            5.68556403e-10 + 6.07307315e-09j,
            3.52904515e-01 + 1.33421124e-03j,
            3.70726944e-01 + 1.59398378e-03j,
            3.91692827e-01 + 1.84823778e-03j,
            4.17420612e-01 + 2.14144753e-03j,
            5.05978088e-01 + 2.03149319e-03j,
            5.60022784e-01 + 2.91701970e-03j,
            5.64743880e-01 + 3.02271626e-03j,
            6.08650092e-01 + 3.38340259e-03j,
            7.27929562e-01 + 2.91272819e-03j,
            7.35090190e-01 + 2.87022736e-03j,
        ]
    )

    evTM_ref = np.array(
        [
            1.13659014e-08 + 2.33360096e-09j,
            3.76427700e-01 + 1.69522962e-03j,
            3.89572593e-01 + 2.15363261e-03j,
            4.17079628e-01 + 2.28176068e-03j,
            4.20818010e-01 + 2.24666939e-03j,
            5.65247098e-01 + 3.11637502e-03j,
            5.66099970e-01 + 3.06711697e-03j,
            5.83461144e-01 + 3.32122025e-03j,
            5.89726166e-01 + 3.18315367e-03j,
            7.45237846e-01 + 3.26633936e-03j,
            7.63801767e-01 + 3.52970555e-03j,
        ]
    )

    npt.assert_almost_equal(evTE, evTE_ref, decimal=2)
    npt.assert_almost_equal(evTM, evTM_ref, decimal=2)

    fem.postpro_eigenvectors()
    pnts = fem.points_kspace(13)
    # return evTE, evTM
