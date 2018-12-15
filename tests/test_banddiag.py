import numpy as np
from pytheas.tools.plottools import *
from pytheas.material import genmat

from pytheas import periodic2D
from pytheas.banddiag.per2D.femmodel import BandsFEM2D
from pytheas.tools import utils
import numpy.testing as npt

pi = np.pi


def model(verbose=False):
    fem = BandsFEM2D()
    # opto-geometric parameters  -------------------------------------------
    mum = 1e-6  #: flt: the scale of the problem (here micrometers)
    fem.dx = 1
    fem.dy = 1
    fem.lambda0search = 10
    fem.neig = 11
    fem.Nix = 101
    fem.parmesh = 10
    fem.type_des = "elements"
    if verbose:
        fem.getdp_verbose = 4
        fem.gmsh_verbose = 4
        fem.python_verbose = 1
    fem.initialize()
    fem.make_mesh()
    return fem


def test_eigpb(verbose=False):
    fem = model(verbose=verbose)
    genmat.np.random.seed(100)
    mat = genmat.MaterialDensity()  # instanciate
    mat.n_x, mat.n_y, mat.n_z = 2 ** 8, 2 ** 8, 1  # sizes
    mat.xsym = True  # symmetric with respect to x?
    mat.p_seed = mat.mat_rand  # fix the pattern random seed
    mat.nb_threshold = 3  # number of materials
    fem.matprop_pattern = [1.4, 2 - 0.02 * 1j, 3 - 0.01j]  # refractive index values
    mat._threshold_val = np.random.permutation(mat.threshold_val)
    mat.pattern = mat.discrete_pattern
    fem.register_pattern(mat.pattern, mat._threshold_val)
    fem.kx, fem.ky = 0, 0
    fem.pola = "TE"
    # fem.initialize()
    fem.compute_solution()
    _ = fem.postpro_eigenvalues()
    evTE = np.sort(_) * fem.dx / (2 * pi)
    fem.pola = "TM"
    fem.update_params()
    fem.compute_solution()
    _ = fem.postpro_eigenvalues()
    evTM = np.sort(_) * fem.dx / (2 * pi)

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
    # return evTE, evTM
