

import numpy as np
from pytheas.tools.plottools import *
from pytheas.material import genmat

from pytheas import periodic2D
from pytheas.periodic2D import FemModel, utils


def model(verbose=False):
    fem = FemModel()
    # opto-geometric parameters  -------------------------------------------
    mum = 1e-6  #: flt: the scale of the problem (here micrometers)
    fem.d = 0.4 * mum  #: flt: period
    fem.h_sup = 1. * mum  #: flt: "thickness" superstrate
    fem.h_sub = 1. * mum  #: flt: "thickness" substrate
    fem.h_layer1 = 0.1 * mum  #: flt: thickness layer 1
    fem.h_layer2 = 0.1 * mum  #: flt: thickness layer 2
    fem.h_des = 0.4 * mum  #: flt: thickness layer design
    fem.h_pmltop = 1. * mum  #: flt: thickness pml top
    fem.h_pmlbot = 1. * mum  #: flt: thickness pml bot
    fem.a_pml = 1  #: flt: PMLs parameter, real part
    fem.b_pml = 1  #: flt: PMLs parameter, imaginary part
    fem.eps_sup = 1  #: flt: permittivity superstrate
    fem.eps_sub = 11  #: flt: permittivity substrate
    fem.eps_layer1 = 1  #: flt: permittivity layer 1
    fem.eps_layer2 = 1  #: flt: permittivity layer 2
    fem.eps_des = 1  #: flt: permittivity layer design
    fem.lambda0 = 0.6 * mum  #: flt: incident wavelength
    fem.theta_deg = 0.  #: flt: incident angle
    fem.pola = "TE"  #: str: polarization (TE or TM)
    fem.lambda_mesh = 0.6 * mum  #: flt: incident wavelength
    #: mesh parameters, correspond to a mesh size of lambda_mesh/(n*parmesh),
    #: where n is the refractive index of the medium
    fem.parmesh_des = 11
    fem.parmesh = 11
    fem.parmesh_pml = fem.parmesh * 2 / 3
    fem.type_des = "elements"
    if verbose:
        fem.getdp_verbose = 4
        fem.gmsh_verbose = 4
        fem.python_verbose = 1
    fem.initialize()
    fem.make_mesh()
    return fem

def aatest_per2D():
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

    fem.register_pattern(mat.pattern, mat._threshold_val)
    fem.compute_solution()
    effs_TE = fem.diffraction_efficiencies()
    fem.pola = "TM"
    fem.compute_solution()
    effs_TM = fem.diffraction_efficiencies()
    print("efficiencies TM", effs_TM)
    rc = 0
    assert rc == 0

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
    fem.open_gmsh_gui()
    rc = 0
    assert rc == 0
