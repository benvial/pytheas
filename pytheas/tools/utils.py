import datetime
import numpy as np
from ..tools import femio
from ..material import genmat


def normalize(x):
    return (x - x.min()) / (x.max() - x.min())


def between_range(x, xmin, xmax):
    return (xmax - xmin) * x + xmin


def generate_ID(self):
    return datetime.now().strftime("%Y_%m_%d_%H_%M_%s_%f")


def refine_mesh(
    fem, mat, lc_des=None, par=None, nmax=5, periodic_x=False, periodic_y=False
):
    mat_tmp = mat
    mat_tmp.ratio_filter = [20, 20, 20]
    pattern = mat_tmp.pattern
    nodes, els, des = fem.get_mesh_info()
    if not lc_des:
        lc_des = 1 * fem.lambda_mesh / (fem.parmesh_des * np.sqrt(fem.eps_des.real))
    # par = [[0.5, 0.4, 0.05], [0.5, 0.4, 0.2], [1, 1, 1]]
    if not par:
        # par = [[0.1], [0.4], [1]]
        par = [[1, 0.8, 0.7], [0.2, 0.15, 0.1], [0.6, 0.8, 1]]
    it = 0
    for smooth_factor, fc_min, fc_max in zip(*par):
        # break
        lc_min, lc_max = lc_des * fc_min, lc_des * fc_max
        # pattern1 = mat.filtered_pattern
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
        # grad_pat_norm = genmat.make_discrete_pattern(grad_pat_norm, np.linspace(0,1,2))
        grad_pat_norm = between_range(grad_pat_norm, lc_min, lc_max)

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
        bgm = "-bgm " + path_size_mesh
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
