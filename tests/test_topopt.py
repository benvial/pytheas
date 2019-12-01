import numpy.testing as npt
from testutils import *

import os
import numpy as np

from pytheas.material import genmat
from pytheas import Scatt2D
import importlib

from pytheas.optim import topopt
from pytheas.optim import TopOpt


# plt.close("all")
pi = np.pi

# np.random.seed(1234)

# ##########################################
# ##############  FEM PARAMETERS  ##########
# ##########################################

#  instanciate fem model ------------------------------------------
fem = Scatt2D()
fem.analysis = "direct"
fem.lambda0 = 1.0
fem.lambda_mesh = fem.lambda0
fem.pola = "TE"
fem.theta_deg = 0.0
fem.h_pml = fem.lambda0
fem.space2pml_L, fem.space2pml_R = fem.lambda0 * 1, fem.lambda0 * 1
fem.space2pml_T, fem.space2pml_B = fem.lambda0 * 1, fem.lambda0 * 1
# material permittivities
tandelta = 0.0
eps_interp = np.array([1, 1.4 ** 2]) * (1 - tandelta * 1j)
fem.eps_interp = eps_interp
fem.eps_des = max(eps_interp)
fem.eps_host = 1.0
fem.parmesh_des = 4
fem.parmesh = 4
fem.parmesh_pml = fem.parmesh * 2 / 3
fem.type_des = "elements"
fem.quad_mesh_flag = True
fem.Nix = 300
fem.Niy = 100
fem.inclusion_flag = False

# size of design box
fem.hx_des = fem.lambda0 * 2
fem.hy_des = fem.lambda0 * 2

# source position
fem.ls_flag = False
fem.xs = 0
fem.ys = +fem.hy_des / 2 + fem.lambda0 / 1

# target
fem.target_flag = True
fem.x_target = 0
fem.y_target = -fem.hy_des / 2 - 0.5 * fem.lambda0
fem.r_target = 0.05

xtar0 = fem.hx_des / 2 * 1.5

####
fem.initialize()
fem.make_mesh()
nvar = len(fem.des[0])

# ##########################################
# #########  OPTIMIZATION PARAMETERS  ######
# ##########################################
to = TopOpt(fem)
to.type_des = fem.type_des
to.algorithm = topopt.nlopt.LD_MMA
to.typeopt = "max"  # type of optimization "min" or "max"
to.pmin = 0  # minimum value
to.pmax = 1  # maximum value
to.m = 1  # interpolation order eps=(eps_min-eps_max)*x^m-eps_min
to.ptol_rel = 1.0e-6
to.ftol_rel = 1.0e-12
to.stopval = None
to.maxeval = 1  # maximum of function evaluation
to.Nitmax = 2  # maximum number of global iterations
to.N0 = 0  # initial global iterations
# to.beta = 1  # projection parameter
lmin = np.min([fem.hx_des, fem.hy_des])
to.rfilt = 0.05 * lmin  # filter radius
to.filt_weight = "gaussian"
to.dg_dp = 0
to.eps_interp = eps_interp
to.log_opt = False

to.plotconv = True
to.force_xsym = False

to.dp = 1e-7
to.m = 1

ratio_hdes = fem.hy_des / fem.hx_des
n_x = 101
n_y = 100  # int(n_x * ratio_hdes) +1
n_z = 1

to.n_x, to.n_y, to.n_z = n_x, n_y, n_z
# genmat.np.random.seed(1234)
mat = genmat.MaterialDensity()
mat.n_x, mat.n_y, mat.n_z = n_x, n_y, n_z
rat = 20
mat.ratio_filter = (rat, rat * ratio_hdes, 1)
mat.nb_threshold = len(eps_interp)
mat.xsym = True
mat.ysym = False
mat.p_seed = mat.mat_rand


def run_fem(p, lambda_range, sens_ana=False, filt=True, proj=True, rmdata=True):

    fem.adjoint = sens_ana
    G, S = [], []
    for fem.lambda0 in lambda_range:
        fem.initialize()
        fem.make_mesh()
        epsilon = to.make_epsilon(p, filt=filt, proj=proj)
        fem.path_pos = fem.make_eps_pos(fem.des[0], epsilon)
        fem.compute_solution()
        goal = to.get_objective()
        if sens_ana:
            sens = to.get_sensitivity(p, filt=filt, proj=proj)
            # fem.open_gmsh_gui()
        else:
            sens = 0
        G.append(goal)
        S.append(sens)
    # if rmdata:
    #     fem.rm_tmp_dir()
    return G, S, fem


# #########################################
# ########  OBJECTIVE FUNCTION       ######
# #########################################


def f_obj(p, grad, filt=True, proj=True, save_history=False):

    # enforce x symetry

    if to.force_xsym:
        p = to.make_xsym(p)

    sens_ana = np.size(grad) > 0
    # nl = 1
    # lambda_range = np.linspace(1, 1.05, nl)
    # grad = np.zeros_like(grad)
    goal, sens = 0, 0
    objpola = []
    x_position = [fem.x_target]  # np.linspace(0, 0.5, n_source)*fem.hx_des/2
    for iobj, fem.x_target in enumerate(x_position):
        lambda_range = [fem.lambda0]
        print("source position: ", fem.x_target)
        G, S, _ = run_fem(p, lambda_range, sens_ana=sens_ana, filt=filt, proj=proj)
        g = np.mean(G, axis=0)
        goal += g
        print(("   objective =  %s " % g))
        objpola.append(g)
        #
        if sens_ana:
            sens += np.mean(S, axis=0)
            if to.force_xsym:
                to.make_xsym(sens)

    grad[:] = sens
    to.obj_history.append(objpola)
    # print(to.obj_history)
    to.param_history.append(p)
    to.tot_obj_history.append(goal)
    to.Nit_tot += 1
    to.Nit_loc += 1

    print("TOTAL")
    print(("   objective =  %s " % goal))
    print("x" * 100)
    return goal


def main_mo(w1, p0=None, save_history=False):
    def f(p, grad, filt=True, proj=True):
        return f_obj(p, grad, filt=filt, proj=proj, save_history=save_history)

    # ##### MAIN OPTIMIZATION LOOP ############
    popt, opt_f, opt = to.main_loop_topopt(f, p0)
    print("optimum at ", popt)
    print("with value  = ", opt_f)
    print(popt)

    if to.threshold_final:
        print("\n")
        print("Final design")
        print("#" * 60)
        if to.force_xsym:
            popt = to.make_xsym(popt)
        popt_filt, _ = to.filter_param(popt, grad=False)
        popt = to.get_threshold_design(popt_filt)
        opt_f = f(popt, np.array([]), filt=False, proj=False)
        print("optimum at ", popt)
        print("with value  = ", opt_f)
    return popt, opt_f, to, fem


def test_topopt():
    # define initial density p0
    random_coef = 1  # np.random.rand()
    p0 = random_coef * to.random_pattern(mat)
    # p0 = 1 * np.random.rand(nvar)
    # w1 = np.linspace(0, 1, 21)
    # P = f_parallel([0], p0=p0, save_history=True)
    out = main_mo(0, p0=p0, save_history=False)
