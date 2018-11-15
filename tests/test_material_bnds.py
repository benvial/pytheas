from pytheas.material.refractiveindex import *


def test_bounds():
    yamlFile = "main/Au/Johnson.yml"
    bounds = getRange(yamlFile)
    assert bounds == (0.1879, 1.937)


def test_vals():
    yamlFile = "main/Au/Johnson.yml"
    ncomplex = get_complex_index(0.5, yamlFile)
    assert np.allclose(ncomplex, 0.97112 - 1.873672j)
