import numpy as np
from pytheas import genmat
from pytheas import utils


def pattern(xsym=False, ysym=False, threeD=False):
    genmat.np.random.seed(100)
    mat = genmat.MaterialDensity()  # instanciate
    mat.n_x, mat.n_y, mat.n_z = 2 ** 8, 2 ** 8, 1  # sizes
    if threeD:
        mat.n_z = mat.n_x
    mat.xsym = xsym  # symmetric with respect to x?
    mat.ysym = ysym  # symmetric with respect to y?
    mat.p_seed = mat.mat_rand  # fix the pattern random seed
    mat.nb_threshold = 3  # number of materials
    mat._threshold_val = np.random.permutation(mat.threshold_val)
    mat.pattern = mat.discrete_pattern
    return mat


def ref_mesh(fem, mat, **kwargs):
    lc_des = fem.lambda0 / (fem.eps_des.real ** 0.5 * fem.parmesh_des)
    par = [[0.5, 0.4, 0.3], [0.4, 0.3, 0.1], [0.8, 0.9, 1]]
    fem = utils.refine_mesh(fem, mat, lc_des=lc_des, par=par, nmax=10, **kwargs)
    return fem
