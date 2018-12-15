from pytheas.material.refractiveindex import *

mat_test = [
    "other/semiconductor alloys/AlN-Al2O3/Hartnett-6.69.yml",
    "other/commercial plastics/Optorez1330/Sultanova.yml",
    "other/doped crystals/Mg-LiTaO3/Moutzouris-e.yml",
    "other/doped crystals/Nb-RbTiOPO4/Carvajal-gamma.yml",
    "other/resists/Microchem 950.yml",
    "other/mixed gases/air/Ciddor.yml",
    "main/Si/Edwards.yml",
    "main/AgBr/Schroter.yml",
    "organic/CH4N2O - urea/Rosker-e.yml",
    "other/clays/montmorillonite/Querry.yml",
    "other/exotic/metamaterials/Valentine.yml",
    "other/heat transfer fluids/Therminol VP-1/Otanicar.yml",
]


def test_all():
    for yamlFile in mat_test:
        bounds = getRange(yamlFile)
        lambdas = np.linspace(bounds[0], bounds[1], 4)
        ncomplex = get_complex_index(lambdas, yamlFile)


def test_bounds():
    yamlFile = "main/Au/Johnson.yml"
    bounds = getRange(yamlFile)
    assert bounds == (0.1879, 1.937)


def test_vals():
    yamlFile = "main/Au/Johnson.yml"
    ncomplex = get_complex_index(0.5, yamlFile)
    assert np.allclose(ncomplex, 0.97112 - 1.873672j)
