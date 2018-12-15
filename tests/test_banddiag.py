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
            7.39746895e-10 + 4.94481158e-09j,
            3.51816155e-01 + 1.32749370e-03j,
            3.71473380e-01 + 1.60101793e-03j,
            3.92461584e-01 + 1.87406751e-03j,
            4.16682124e-01 + 2.19776755e-03j,
            5.10246232e-01 + 2.07844081e-03j,
            5.62985322e-01 + 3.04380085e-03j,
            5.64954110e-01 + 3.07417776e-03j,
            5.99827366e-01 + 3.36526859e-03j,
            7.25296682e-01 + 2.85387396e-03j,
            7.36953867e-01 + 2.90558327e-03j,
        ]
    )

    evTM_ref = np.array(
        [
            1.18847963e-08 - 7.98118949e-10j,
            3.76843862e-01 + 1.70492358e-03j,
            3.92827392e-01 + 2.23200466e-03j,
            4.15643814e-01 + 2.33824362e-03j,
            4.17611545e-01 + 2.30297337e-03j,
            5.60156432e-01 + 3.06420619e-03j,
            5.66434068e-01 + 3.01957339e-03j,
            5.83496928e-01 + 3.46354189e-03j,
            5.85668625e-01 + 3.29038632e-03j,
            7.48985237e-01 + 3.34867583e-03j,
            7.57608348e-01 + 3.34302387e-03j,
        ]
    )

    npt.assert_almost_equal(evTE, evTE_ref, decimal=3)
    npt.assert_almost_equal(evTM, evTM_ref, decimal=3)
    # return evTE, evTM
