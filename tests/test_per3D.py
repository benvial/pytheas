import numpy as np
from pytheas import Periodic3D
import numpy.testing as npt
from testutils import *


def model(verbose=False):
    fem = Periodic3D()
    fem.rm_tmp_dir()
    fem.lambda0 = 1.1
    fem.theta_deg = 0
    fem.phi_deg = 0
    fem.psi_deg = 0
    fem.period_x = 1.0  #: flt: periodicity in x-direction
    fem.period_y = 1.0  #: flt: periodicity in y-direction
    fem.thick_L1 = 1
    fem.thick_L2 = 0.1
    fem.thick_L3 = 0.05
    fem.thick_L4 = 0.1
    fem.thick_L5 = 0.1
    fem.thick_L6 = 1
    fem.PML_top = 1.4  #: flt: thickness pml top
    fem.PML_bot = 1.4  #: flt: thickness pml bot
    fem.eps_L2 = 1 - 0.0 * 1j  #: flt: permittivity layer 2
    fem.eps_L3 = 4 - 0.2 * 1j  #: flt: permittivity layer 3
    fem.eps_L4 = 1 - 0.0 * 1j  #: flt: permittivity layer 4
    fem.eps_L5 = 1 - 0.0 * 1j  #: flt: permittivity layer 5
    fem.parmesh_des = 6
    fem.parmesh = 5
    fem.parmesh_pml = fem.parmesh * 2 / 3
    fem.N_d_order = 0
    fem.el_order = 2  #: int: order of basis function (1 or 2)
    fem.matprop_pattern = [1.4, 2 - 0.02 * 1j, 3 - 0.01j]  # refractive index values
    if verbose:
        fem.getdp_verbose = 4
        fem.gmsh_verbose = 4
        fem.python_verbose = 1
    fem.initialize()
    fem.make_mesh()

    return fem


def test_per3D(verbose=True):
    fem = model(verbose=verbose)
    mat = pattern(xsym=True, ysym=True)
    fem.register_pattern(mat.pattern, mat._threshold_val)
    fem.compute_solution()
    effs = fem.diffraction_efficiencies()
    print("effs = ", effs)
    effs_ref = {
        "T": np.array([[0.796112]]),
        "R": np.array([[0.18909322]]),
        "Q": 0.01273080344897011,
        "B": 0.9979360252786229,
    }

    for a, b in zip(effs.values(), effs_ref.values()):
        npt.assert_almost_equal(a, b, decimal=2)
