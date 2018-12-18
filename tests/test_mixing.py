from pytheas import mixing
import numpy.testing as npt


def test_mixing():
    incl_vol_frac, eps_incl, eps_host = 0.5, 9, 3
    npt.assert_almost_equal(
        mixing.arithmetic_mean(incl_vol_frac, eps_incl, eps_host), 6, decimal=3
    )
    npt.assert_almost_equal(
        mixing.harmonic_mean(incl_vol_frac, eps_incl, eps_host), 4.5, decimal=3
    )
    npt.assert_almost_equal(
        mixing.maxwell_garnett(incl_vol_frac, eps_incl, eps_host), 5.25, decimal=3
    )
