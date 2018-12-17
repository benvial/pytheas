import numpy as np
from pytheas import genmat
from pytheas import utils


def pattern():
    genmat.np.random.seed(100)
    mat = genmat.MaterialDensity()  # instanciate
    mat.n_x, mat.n_y, mat.n_z = 2 ** 8, 2 ** 8, 1  # sizes
    mat.xsym = True  # symmetric with respect to x?
    mat.p_seed = mat.mat_rand  # fix the pattern random seed
    mat.nb_threshold = 3  # number of materials
    mat._threshold_val = np.random.permutation(mat.threshold_val)
    mat.pattern = mat.discrete_pattern
    return mat


def ref_mesh(fem, mat):
    lc_des = fem.lambda0 / (fem.eps_des.real ** 0.5 * fem.parmesh_des)
    par = [[0.5, 0.4, 0.3], [0.4, 0.3, 0.1], [0.8, 0.9, 1]]
    fem = utils.refine_mesh(fem, mat, lc_des=lc_des, par=par, periodic_x=False, nmax=10)
    return fem
