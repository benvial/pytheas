"""Shared utility functions used in pytheas."""


import numpy as np
from ..tools import femio
from ..material import genmat


def normalize(x):
    """Normalize an array between 0 and 1

        Parameters
        ----------
        x : array-like
            the quantity to be normalized

        Returns
        -------
        x_norm : array-like
            normalized array


    """
    return (x - x.min()) / (x.max() - x.min())


def between_range(x, xmin, xmax):
    return (xmax - xmin) * x + xmin


def refine_mesh(
    fem,
    mat,
    lc_des=None,
    par=None,
    nmax=5,
    periodic_x=False,
    periodic_y=False,
    lc_bnd=None,
):
    mat_tmp = mat
    mat_tmp.ratio_filter = [20, 20, 20]
    pattern = mat_tmp.pattern
    nodes, els, des = fem.get_mesh_info()
    if not lc_des:
        lc_des = 1 * fem.lambda_mesh / (fem.parmesh_des * np.sqrt(fem.eps_des.real))
    if not par:
        par = [[1, 0.8, 0.7], [0.2, 0.15, 0.1], [0.6, 0.8, 1]]
    it = 0
    for smooth_factor, fc_min, fc_max in zip(*par):
        lc_min, lc_max = lc_des * fc_min, lc_des * fc_max
        if fem.dim is 2:
            grad_pat = np.gradient(pattern[:, :, 0])
            grad_pat_norm = np.sqrt(grad_pat[0] ** 2 + grad_pat[1] ** 2).reshape(
                pattern.shape
            )
        else:
            grad_pat = np.gradient(pattern)
            grad_pat_norm = np.sqrt(
                grad_pat[0] ** 2 + grad_pat[1] ** 2 + grad_pat[2] ** 2
            )

        grad_pat_norm = genmat.filter_pattern(
            grad_pat_norm, smooth_factor * mat_tmp.sigma
        )
        grad_pat_norm = normalize(grad_pat_norm)
        thres_grad = 1e-1
        grad_pat_norm[grad_pat_norm > thres_grad] = 1
        grad_pat_norm[grad_pat_norm <= thres_grad] = 0
        grad_pat_norm = genmat.filter_pattern(
            grad_pat_norm, smooth_factor * mat_tmp.sigma
        )

        grad_pat_norm = 1 - normalize(grad_pat_norm)
        grad_pat_norm = between_range(grad_pat_norm, lc_min, lc_max)

        if not lc_bnd:
            lc_bnd = lc_min  # (lc_max + lc_min) / 2
        if periodic_y:
            grad_pat_norm[:, :nmax] = lc_bnd
            grad_pat_norm[:, -nmax:] = lc_bnd
        if periodic_x:
            grad_pat_norm[:nmax, :] = lc_bnd
            grad_pat_norm[-nmax:, :] = lc_bnd
        fdens_grad = fem.make_fdens(grad_pat_norm)
        density_grad = fdens_grad(nodes[1])
        fem.type_des = "nodes"
        path_size_mesh = fem.make_pos(nodes[0], density_grad, "size_mesh")
        fem.type_des = "elements"
        # if it == 0:
        if fem.dim is 3:
            dim = [1, 2, 3]
        else:
            dim = [1, 2]
        femio.mesh_model(
            fem.path_mesh, fem.path_bg_mesh, verbose=fem.gmsh_verbose, dim=dim
        )

        bgm = " -bgm {}".format(path_size_mesh)
        bgm += " -merge {}".format(fem.path_bg_mesh)
        femio.mesh_model(
            fem.path_mesh,
            fem.path_geo,
            other_option=bgm,
            verbose=fem.gmsh_verbose,
            dim=dim,
        )
        nodes, els, des = fem.get_mesh_info()
        # fem.content_mesh = fem.make_mesh_pos(els, nodes)
        fem.content_mesh = femio.get_content(fem.path_mesh)
        it += 1
        if fem.gmsh_verbose:
            print("mesh refinement iteration :", it)
        fem.nodes, fem.els, fem.des = nodes, els, des
    return fem
