import numpy as np
from pytheas.tools.plottools import *
from pytheas.material import genmat

from pytheas import periodic2D
from pytheas.periodic2D import FemModel
from pytheas.tools import utils
import numpy.testing as npt


def model(verbose=False):
    fem = FemModel()
    # opto-geometric parameters  -------------------------------------------
    mum = 1e-6  #: flt: the scale of the problem (here micrometers)
    fem.d = 0.3 * mum  #: flt: period
    fem.h_sup = 1.0 * mum  #: flt: "thickness" superstrate
    fem.h_sub = 1.0 * mum  #: flt: "thickness" substrate
    fem.h_layer1 = 0.1 * mum  #: flt: thickness layer 1
    fem.h_layer2 = 0.1 * mum  #: flt: thickness layer 2
    fem.h_des = 0.4 * mum  #: flt: thickness layer design
    fem.h_pmltop = 0.6 * mum  #: flt: thickness pml top
    fem.h_pmlbot = 0.6 * mum  #: flt: thickness pml bot
    fem.a_pml = 1  #: flt: PMLs parameter, real part
    fem.b_pml = 1  #: flt: PMLs parameter, imaginary part
    fem.eps_sup = 1  #: flt: permittivity superstrate
    fem.eps_sub = 3  #: flt: permittivity substrate
    fem.eps_layer1 = 3  #: flt: permittivity layer 1
    fem.eps_layer2 = 1  #: flt: permittivity layer 2
    fem.eps_des = 9  #: flt: permittivity layer design
    fem.lambda0 = 0.6 * mum  #: flt: incident wavelength
    fem.theta_deg = 0.0  #: flt: incident angle
    fem.pola = "TE"  #: str: polarization (TE or TM)
    fem.lambda_mesh = 0.6 * mum  #: flt: incident wavelength
    fem.quad_mesh_flag = False
    #: mesh parameters, correspond to a mesh size of lambda_mesh/(n*parmesh),
    #: where n is the refractive index of the medium
    fem.parmesh_des = 5
    fem.parmesh = 5
    fem.parmesh_pml = fem.parmesh * 2 / 3
    fem.type_des = "elements"
    if verbose:
        fem.getdp_verbose = 4
        fem.gmsh_verbose = 4
        fem.python_verbose = 1
    fem.initialize()
    fem.make_mesh()

    return fem


def test_per2D(verbose=False):
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

    lc_des = fem.lambda0 / (fem.eps_des.real ** 0.5 * fem.parmesh_des)
    # par = [[0.5, 0.4, 0.15], [0.4, 0.3, 0.1], [0.1, 0.5, 1]]
    par = [[0.5, 0.4, 0.3], [0.4, 0.3, 0.1], [0.8, 0.9, 1]]
    fem = utils.refine_mesh(fem, mat, lc_des=lc_des, par=par, periodic_x=False, nmax=10)

    fem.register_pattern(mat.pattern, mat._threshold_val)
    fem.compute_solution()
    effs_TE = fem.diffraction_efficiencies()
    print("effs_TE = ", effs_TE)
    fem.pola = "TM"
    fem.compute_solution()
    effs_TM = fem.diffraction_efficiencies()
    print("effs_TM = ", effs_TM)

    effs_TE_ref = {
        "R": 0.3528373505073045,
        "T": 0.3274087995895207,
        "Q": 0.3262777713934278,
        "B": 1.006523921490253,
    }
    effs_TM_ref = {
        "R": 0.2966385564434248,
        "T": 0.5843153749078397,
        "Q": 0.1262830646121085,
        "B": 1.007236995963373,
    }

    for a, b in zip(effs_TE.values(), effs_TE_ref.values()):
        npt.assert_almost_equal(a, b, decimal=3)
    for a, b in zip(effs_TM.values(), effs_TM_ref.values()):
        npt.assert_almost_equal(a, b, decimal=3)


def test_ref_mesh():
    fem = model()
    genmat.np.random.seed(100)
    mat = genmat.MaterialDensity()  # instanciate
    mat.n_x, mat.n_y, mat.n_z = 2 ** 7, 2 ** 7, 1  # sizes
    mat.xsym = True  # symmetric with respect to x?
    mat.p_seed = mat.mat_rand  # fix the pattern random seed
    mat.nb_threshold = 3  # number of materials
    fem.matprop_pattern = [1.4, 4 - 0.02 * 1j, 2]  # refractive index values
    mat._threshold_val = np.random.permutation(mat.threshold_val)
    mat.pattern = mat.discrete_pattern
    fem = utils.refine_mesh(fem, mat)
    fem.register_pattern(mat.pattern, mat._threshold_val)
    fem.compute_solution()
    rc = 0
    assert rc == 0
