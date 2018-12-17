import numpy as np
from pytheas import HighContrast2D
import numpy.testing as npt
from testutils import *


def square(a, x0, y0):
    nt = 360 * 4
    theta = np.linspace(0, 2 * np.pi, nt)
    eps = 1e-6
    rt = [
        np.min(
            [0.5 * a / (eps + np.abs(np.cos(t))), 0.5 * a / (eps + np.abs(np.sin(t)))]
        )
        for t in theta
    ]
    rt = np.array(rt)
    x = rt * np.cos(theta)
    y = rt * np.sin(theta)
    points = x + x0, y + y0
    return points


def model(verbose=False):
    fem = HighContrast2D()
    fem.rm_tmp_dir()
    fem.eps_incl = 50 - 1j
    fem.eps_host = 2
    fem.parmesh = 10
    fem.inclusion_flag = True
    fem.inclusion_filename_ = "inclusion.geo"
    fem.neig = 20
    if verbose:
        fem.getdp_verbose = 4
        fem.gmsh_verbose = 4
        fem.python_verbose = 1
    fem.initialize()
    fem.make_inclusion(square(0.8, 0, 0))
    fem.make_mesh()
    return fem


def test_hom(verbose=False):
    fem = model(verbose=verbose)
    fem.compute_modes()
    fem.compute_solution()
    eps_eff = fem.compute_epsilon_eff()
    fem.compute_modes()
    mu_eff = fem.postpro_effective_permeability()

    kref = np.array([0.7, 0.9])
    mu_eff_ref = np.array([2.63250092 - 0.15751293j, -0.71608055 - 0.11285409j])
    npt.assert_almost_equal(mu_eff(kref), mu_eff_ref, decimal=3)
    return fem
