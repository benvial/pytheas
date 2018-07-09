
from pytheas.tools import femio
import numpy as np
from pytheas.tools.plottools import *
from pytheas.material import genmat


pi = np.pi

# # refine mesh ------------------------------------------


def normalize(x):
    return (x - x.min()) / (x.max() - x.min())


def between_range(x, xmin, xmax):
    return (xmax - xmin) * x + xmin


def refine_mesh(fem, mat, lc_des=None, par=None, periodic=False):
    mat_tmp = mat
    mat_tmp.ratio_filter = [20, 20, 20]
    pattern = mat_tmp.pattern
    nodes, els, des = get_mesh_info(fem)
    if not lc_des:
        lc_des = 1 * fem.lambda_mesh / \
            (fem.parmesh_des * np.sqrt(fem.eps_des.real))
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
            grad_pat_norm = np.sqrt(
                grad_pat[0] ** 2 + grad_pat[1] ** 2).reshape(pattern.shape)
        else:
            grad_pat = np.gradient(pattern)
            grad_pat_norm = np.sqrt(
                grad_pat[0] ** 2 + grad_pat[1] ** 2 + grad_pat[2] ** 2)

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
        # for i in range(np.shape(grad_pat)[-1]): plt.clf, plt.imshow(grad_pat_norm[:,:,i]), plt.pause(0.01)

        if periodic:
            lc_bnd = (lc_max + lc_min) / 2
            grad_pat_norm[:, 0] = lc_bnd
            grad_pat_norm[:, -1] = lc_bnd
            grad_pat_norm[0, :] = lc_bnd
            grad_pat_norm[-1, :] = lc_bnd

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
        femio.mesh_model(fem.path_mesh, fem.path_bg_mesh,
                         verbose=fem.gmsh_verbose, dim=dim)
        bgm = "-bgm " + path_size_mesh
        femio.mesh_model(
            fem.path_mesh, fem.path_geo, other_option=bgm, verbose=fem.gmsh_verbose,  dim=dim)
        nodes, els, des = get_mesh_info(fem)
        fem.content_mesh = fem.make_mesh_pos(els, nodes)
        # plt.clf()
        # plt.imshow(np.fliplr(grad_pat_norm[:, :, 0]).T)
        # plt.colorbar()
        # plt.pause(0.1)
        # fem.open_gmsh_gui(pos_list=["size_mesh.pos"])
        it += 1
        if fem.gmsh_verbose:
            print("mesh refinement iteration :", it)
        fem.nodes, fem.els, fem.des = nodes, els, des
    return fem


def get_mesh_info(fem):
    # get nodes and elements and their IDs in the design domain
    nodes = fem.get_design_nodes()
    els = fem.get_design_elements()
    nodes_ID, nodes_coords = nodes
    els_ID, els_coords, els_nodes_ID, geom_ID_dom = els
    xnodes, ynodes, znodes = nodes_coords.T
    xels, yels, zels = els_coords.T
    if fem.type_des is "elements":
        des_ID, des_coords = els_ID, els_coords
        des = els_ID, els_coords
    elif fem.type_des is "nodes":
        des_ID, des_coords = nodes_ID, nodes_coords
        des = nodes_ID, nodes_coords
    return nodes, els, des


def transform_eff(eff):
    effa = np.array(eff).T
    R = [e["R"] for e in eff]
    T = [e["T"] for e in eff]
    Q = [e["Q"] for e in eff]
    B = [e["B"] for e in eff]
    return R, T, Q, B


def normalize_components(a):
    an = np.zeros_like(a)
    for i, c in enumerate(a):
        an[i, :] = (c - c.min()) / (c.max() - c.min())
    return an


def plot_spectrum(lambda_range, eff):
    R, T, Q, B = transform_eff(eff)
    xplot = lambda_range
    # plt.plot(xplot, np.sum(R, axis=1), label='R', color='#0077b3')
    # plt.plot(xplot, np.sum(T, axis=1), label='T', color='#cc3300')
    plt.plot(xplot, R, label="R", color="#0077b3")
    plt.plot(xplot, T, label="T", color="#cc3300")
    plt.plot(xplot, Q, label="A", color="#77b300")
    plt.plot(xplot, B, label="B", color="#666666")
    plt.title("diffraction efficiencies spectrum")
    plt.legend(loc=0)
    plt.xlabel(r"Wavelength $\lambda_0$ (nm)")
    plt.ylabel(r"Diffraction efficiencies")


def meromorphic_expansion(
    fem, eigval, deg_poly, Nlambda, lambdai_min=None, lambdai_max=None
):
    # ###------ meromorphic expansion of r and t--------
    fem.cplx_effs = True
    reson = eigval / (2 * pi)
    lamb_res = 1 / reson.real
    Q_res = 0.5 * np.abs(reson.real / reson.imag)
    # reson = reson[(lamb_res > 1) & (lamb_res < 2)]
    n_res = len(reson)
    Ni = deg_poly + 1 + n_res
    # lamb_poly = []
    # for i in range(n_res - 1):
    #     lamb_poly.append((lamb_res[i] + lamb_res[i + 1]) / 2)
    # lamb_poly = np.random.choice(lamb_poly, deg_poly + 1)
    if not lambdai_min:
        lambdai_min = lamb_res.min() * 0.9
    if not lambdai_max:
        lambdai_max = lamb_res.max() * 1.1
    # lamb_poly = np.random.uniform(lambdai_min, lambdai_max, deg_poly + 1)
    # lamb_poly = np.linspace(lambdai_min, lambdai_max, deg_poly + 1)
    lamb_poly = 1 / np.linspace(1 / lambdai_min, 1 / lambdai_max, deg_poly + 1)
    lambdai = np.hstack((lamb_poly, lamb_res))
    # lambdai = np.linspace(lambdai_min, lambdai_max, Ni)
    f = []
    # oi = 2*pi*fem.cel/lambdai
    oi = 1 / lambdai
    for o in reson:
        f.append(1 / (oi - o))
    for n in range(deg_poly + 1):
        f.append(oi ** n)

    param = lambdai
    eff = lambda_sweep(fem, param)
    # FIXME: make the expansion for every order (this does NOT work for more than one propagating order)
    Rparam = [np.sum(e["R"]) for e in eff]
    Tparam = [np.sum(e["T"]) for e in eff]
    inv_f = np.linalg.inv(np.array((f)).T)
    coefs_R = np.dot(inv_f, (Rparam))
    coefs_T = np.dot(inv_f, (Tparam))
    lambda_mod = np.linspace(lambdai_min, lambdai_max, Nlambda)
    o_mod = 1 / lambda_mod
    R_mod, T_mod = 0, 0
    i = 0
    for o in reson:
        R_mod += coefs_R[i] / (o_mod - o)
        T_mod += coefs_T[i] / (o_mod - o)
        i += 1
    for n in range(deg_poly + 1):
        R_mod += coefs_R[i] * o_mod ** n
        T_mod += coefs_T[i] * o_mod ** n
        i += 1
    Rparam = np.array(Rparam).ravel()
    Tparam = np.array(Tparam).ravel()
    return param, Rparam, Tparam, lambda_mod, R_mod, T_mod


def lambda_sweep(fem, param):
    fem.analysis = "diffraction"
    eff = []
    for par in param:
        fem.lambda0 = par
        fem.initialize()
        fem.compute_solution()
        eff.append(fem.diffraction_efficiencies())
    return eff


def plot_epsilon(fem, pattern, matprop):
    cmap_cont = sns.cubehelix_palette(8, start=0.5, rot=-0.75, as_cmap=True)
    cmap = cmap_discretize(cmap_cont, nmat)
    extent = [-fem.d / 2, fem.d / 2, fem.h_layer1, fem.h_layer1 + fem.h_des]
    plt.imshow(
        np.flipud(pattern[:, :, 0].T),
        interpolation=None,
        cmap=cmap,
        vmin=0,
        vmax=1,
        extent=extent,
    )
    cbar = plt.colorbar(ticks=np.linspace(0, 1, nmat))
    plt.clim(-1 / (nmat + 1), 1 + 1 / (nmat + 1))
    # cbarlabels = ['%01f %+01fi' % (i.real, i.imag) for i in matprop]
    # cbarlabels =  [r'%.2f %+.2f $10^{-3} i$' % (i.real, 1e3*i.imag) for i in matprop]
    # cbarlabels = [r'%.3f %+.3fi' % (i.real, i.imag) for i in matprop]
    cbarlabels = [
        r"%.1f %+.1f $10^{-2} i$" % ((i ** 2).real, 1e2 * (i ** 2).imag)
        for i in matprop
    ]
    cbar.ax.set_yticklabels(cbarlabels)
    plt.xlabel("d = {0:.4}".format(fem.d))
    plt.ylabel("h = {0:.4}".format(fem.h_des))
