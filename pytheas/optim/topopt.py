import numpy as np
import nlopt
import numpy as np
import scipy as sc
from pytheas.tools.plottools import *

# from scipy.interpolate import splev, splrep
from scipy.interpolate import PchipInterpolator  # mono_cubic_interp


def normalize(x0):
    if x0.min() == x0.max():
        return x0
    else:
        return (x0 - x0.min()) / (x0.max() - x0.min())


def norm_vec(x, y):
    return np.sqrt((x) ** 2 + (y) ** 2)


def simp_(p, val, m=1):
    nthres = len(val)
    # k = min(self.nthres - 1, 3)
    if nthres == 2:
        out = (val[1] - val[0]) * p ** m + val[0]
    else:
        p_interp = np.linspace(0, 1, nthres)
        # tsimp_re = splrep(self.p_interp, self.eps_interp.real, k=k)
        # tsimp_im = splrep(self.p_interp, self.eps_interp.imag, k=k)
        # out = splev(p, tsimp_re) + 1j*splev(p, tsimp_im)
        tsimp_re = PchipInterpolator(p_interp, val.real)
        tsimp_im = PchipInterpolator(p_interp, val.imag)
        out = tsimp_re(p ** m) + 1j * tsimp_im(p ** m)
    return out


class TopologyOptimization:
    """A class for topology optimization

    """

    def __init__(self, fem):
        self.obj_history = []
        self.param_history = []
        self.tot_obj_history = []
        self.fem = fem

    ###########################################
    ##########  OPTIMIZATION PARAMETERS  ######
    ###########################################

    algorithm = nlopt.LD_MMA
    type_des = "elements"
    typeopt = "min"  # type of optimization "min" or "max"
    pmin = 0.0  # minimum value
    pmax = 1.0  # maximum value
    m = 1.0  # interpolation order eps=(eps_min-eps_max)*x^m-eps_min
    ptol_rel = 1.0e-8
    ftol_rel = 1.0e-16
    stopval = None
    maxeval = 50  # maximum of function evaluation
    Nitmax = 8  # maximum number of global iterations
    Nit = 0  # initialize global iteration number
    Nit_loc = 0
    N0 = 0
    Nit_tot = 0
    beta = 1  # projection parameter
    rfilt = 0  # filter radius
    filt_weight = "gaussian"
    plotconv = False
    dg_dp = 0
    eps_interp = np.array([1, 5])
    threshold_final = True

    log_opt = False
    dp = 1e-11
    n_x, n_y, n_z = 100, 100, 1
    force_xsym = False
    # obj_history = []
    # param_history = []

    @property
    def nthres(self):
        return len(self.eps_interp)

    @property
    def thres(self):
        return np.linspace(0, 1, self.nthres + 1)[1 : (self.nthres)]

    @property
    def p_interp(self):
        return np.linspace(0, 1, self.nthres)

    @property
    def grid(self):
        xdes, ydes, zdes = self.fem.des[1].T
        x1 = np.linspace(min(xdes), max(xdes), self.n_x)
        y1 = np.linspace(min(ydes), max(ydes), self.n_y)
        return x1, y1

    def simp(self, p):
        return simp_(p, self.eps_interp, m=self.m)

    def proj(self, x):
        p = 0
        for thres in self.thres:
            ptmp = np.tanh(self.beta * (x - thres)) / (np.tanh(thres * self.beta))
            # p += normalize(ptmp)
            p += ptmp
        p = (p + self.nthres - 1) / (2 * self.nthres - 2)

        return p
        # return x

    def grad_simp(self, p):
        f = self.simp(p)
        df = self.simp(p + self.dp)
        return (df - f) / self.dp

    def grad_proj(self, p):
        f = self.proj(p)
        df = self.proj(p + self.dp)
        return (df - f) / self.dp

    # Filter
    def neighbourhood(self, xe, ye, xn, yn):
        return norm_vec(xn - xe, yn - ye) <= self.rfilt

    def weight_filt(self, x, y, xe, ye):
        D = norm_vec(x - xe, y - ye)
        if self.filt_weight is "gaussian":
            w = np.exp(-0.5 * (2 * D / self.rfilt) ** 2)
        elif self.filt_weight is "linear":
            w = self.rfilt - D
        elif self.filt_weight is "constant":
            w = 1.0
        else:
            raise ValueError(
                "Wrong weight type '{0}'. Available weights are 'gaussian', 'linear' or 'constant'".format(
                    self.filt_weight
                )
            )
        return w

    def filter_param(self, p):
        xdes, ydes, zdes = self.fem.des[1].T
        if self.rfilt == 0:
            pfilt = p
        else:
            pfilt = []
            for xe, ye in zip(xdes, ydes):
                Ne = self.neighbourhood(xe, ye, xdes, ydes)
                xi, yi, pi = xdes[Ne], ydes[Ne], p[Ne]
                wi = self.weight_filt(xi, yi, xe, ye)
                ptmp = sum(wi * pi) / sum(wi)
                pfilt.append(ptmp)
        return np.array(pfilt)
        # return np.array(p)

    def grad_filter_param(self, p):
        xdes, ydes, zdes = self.fem.des[1].T
        f = self.filter_param(p)
        df = self.filter_param(p + self.dp)
        return (df - f) / self.dp

    def make_xsym(self, qt, interp_method="cubic"):
        qt = self.mesh2grid(qt, interp_method=interp_method)
        qt = (qt + np.fliplr(qt)) / 2
        return self.grid2mesh(qt, interp_method=interp_method)

    def random_pattern(self, mat):
        im = mat.filtered_pattern[:, :, 0].T
        x_grid = np.linspace(self.grid[0][0], self.grid[0][-1], mat.n_x)
        y_grid = np.linspace(self.grid[1][0], self.grid[1][-1], mat.n_y)
        x0 = self.grid2mesh(im, grid=(x_grid, y_grid))
        return normalize(x0)

    def grid2mesh(self, val_grid, grid=None, interp_method="cubic"):
        xdes, ydes, zdes = self.fem.des[1].T
        if grid:
            x_grid, y_grid = grid
        else:
            x_grid, y_grid = self.grid
        f = sc.interpolate.interp2d(x_grid, y_grid, val_grid, kind=interp_method)
        x0 = [f(x, y) for x, y in zip(xdes, ydes)]
        return np.array(x0).ravel()

    def mesh2grid(self, valdes, interp_method="cubic"):
        xdes, ydes, zdes = self.fem.des[1].T
        points = np.vstack((xdes, ydes)).T
        x_grid, y_grid = self.grid
        valdes.flags.writeable = True
        xg, yg = np.meshgrid(x_grid, y_grid)
        v = sc.interpolate.griddata(points, valdes, (xg, yg), method=interp_method)
        return v

    #

    def make_epsilon(self, p, filt=True, proj=True):
        xdes, ydes, zdes = self.fem.des[1].T
        if filt:
            p_filt = self.filter_param(p)
        else:
            p_filt = p
        if proj:
            p_proj = self.proj(p_filt)
        else:
            p_proj = p_filt
        epsilon = self.simp(p_proj)
        return epsilon

    def depsilon_dp(self, p, filt=True, proj=True):
        if filt:
            p_filt = self.filter_param(p)
            dpfilt_dp = self.grad_filter_param(p)
        else:
            p_filt = p
            dpfilt_dp = np.ones_like(p)
        if proj:
            p_proj = self.proj(p_filt)
            dpproj_dpfilt = self.grad_proj(p_filt)
        else:
            p_proj = p_filt
            dpproj_dpfilt = np.ones_like(p)
        depsilon_dpproj = self.grad_simp(p_proj)
        return depsilon_dpproj * dpproj_dpfilt * dpfilt_dp

    def get_objective(self):
        goal = self.fem.get_objective()
        # if self.log_opt:
        #     goal = np.log(goal)
        return goal

    def get_adjoint(self):
        return self.fem.get_adjoint()

    def get_deq_deps(self, interp_method="cubic"):
        deq_deps = self.fem.get_deq_deps()

        if self.fem.pola is "TM":
            x_grid, y_grid = self.grid

            deq_deps_x, deq_deps_y = deq_deps
            deq_deps_x = self.mesh2grid(deq_deps_x, interp_method=interp_method)
            deq_deps_x_x = np.gradient(deq_deps_x.T)[0] / np.gradient(x_grid)[0]

            deq_deps_y = self.mesh2grid(deq_deps_y, interp_method=interp_method)
            deq_deps_y_y = np.gradient(deq_deps_y.T)[1] / np.gradient(y_grid)[0]

            deq_deps = deq_deps_x_x.T + deq_deps_y_y.T

            deq_deps_re = self.grid2mesh(deq_deps.real)
            deq_deps_im = self.grid2mesh(deq_deps.imag)
            deq_deps = deq_deps_re + 1j * deq_deps_im
        return deq_deps

    def get_sensitivity(self, p, filt=True, proj=True, interp_method="cubic"):
        adjoint = self.get_adjoint()
        deq_deps = self.get_deq_deps(interp_method=interp_method)
        sens = self.dg_dp + 1 * np.real(
            adjoint * deq_deps * self.depsilon_dp(p, filt=filt, proj=proj)
        )
        # if self.log_opt:
        #     p_min = 1e-3
        #     p_thres = np.copy(p)
        #     p_thres[p <= p_min] = p_min
        #     sens /= p_thres
        return sens

    def plot_design(
        self, ax, varplot, typeplot="interp", cmap=None, extent=None, **kwargs
    ):
        if typeplot is "tri":
            xdes, ydes, zdes = self.fem.des[1].T
            triang = matplotlib.tri.Triangulation(xdes, ydes)
            xmid = xdes[triang.triangles].mean(axis=1)
            ymid = ydes[triang.triangles].mean(axis=1)
            cf = ax.tripcolor(triang, varplot, shading="flat", cmap=cmap)
        elif typeplot is "interp":
            varplot = self.mesh2grid(varplot, **kwargs)
            cf = ax.imshow(np.flipud(varplot), cmap=cmap, extent=extent)
        else:
            raise TypeError("Wrong typeplot specified: choose between interp and tri")
        cbar = plt.colorbar(cf, fraction=0.046, pad=0.04)
        # ax.set_title(r'Permittivity $\varepsilon$ in the design domain')
        ax.axis("image")
        ax.axis("off")

    def plot_convergence(self, ax):
        if self.Nit_tot > 20:
            styleplot = "-"
        else:
            styleplot = "-o"

        x_obj = list(range(self.Nit_tot))
        obj_array = np.array(self.obj_history).T
        for obj in obj_array:
            ax.plot(x_obj, obj, styleplot)
            ax.fill_between(x_obj, obj, 0.0, alpha=0.2)
        ax.plot(x_obj, self.tot_obj_history, "--", color="gray")
        ax.set_title(
            "Convergence: "
            + "global iteration = "
            + str(self.Nit)
            + ", current iteration = "
            + str(self.Nit_tot)
        )
        ax.set_xlabel(r"iteration number $N$")
        ax.set_ylabel(r"objective $\varphi$")
        # ax4.grid(True)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.axis("tight")

    def plot_while_solving(self, varplot, title="", **kwargs):
        # print("beta = ", beta)
        plt.clf()
        ax1 = plt.subplot(211, aspect="equal")
        self.plot_design(ax1, varplot, **kwargs)
        ax1.set_title(title)
        ax2 = plt.subplot(212)
        self.plot_convergence(ax2)
        if self.Nit_tot > 1:
            ax2.set_xlim((0, self.Nit_tot - 1))
        # ax2.set_ylim((0,1))
        plt.tight_layout()
        # plt.draw()
        plt.show()
        plt.pause(0.01)

    def main_loop_topopt(self, f_obj, p0):
        nvar = len(p0)
        print("\n")
        print("#" * 60)
        print("Topology optimization with %s variables" % nvar)
        print("#" * 60)
        # instanciation opt object, optimization algorithm
        opt = nlopt.opt(self.algorithm, nvar)
        # opt.verbose=1
        opt.set_lower_bounds(self.pmin)
        opt.set_upper_bounds(self.pmax)
        # opt.set_min_objective(f_obj)
        if self.typeopt is "max":
            opt.set_max_objective(f_obj)
        elif self.typeopt is "min":
            opt.set_min_objective(f_obj)
        else:
            raise TypeError("Wrong typeopt specified: choose between max and min")
        # opt.add_inequality_constraint(lambda x,grad: myconstraint(x,grad,2,0), 1e-8)
        # opt.add_inequality_constraint(lambda x,grad: myconstraint(x,grad,-1,1), 1e-8)
        opt.set_maxeval(self.maxeval)
        if self.stopval:
            opt.set_stopval(self.stopval)
        opt.set_xtol_rel(self.ptol_rel)
        opt.set_ftol_rel(self.ftol_rel)
        self.opt = opt
        ################################################################
        ############### OPTIMIZATION - global iterations ###############
        ################################################################
        for self.Nit in range(self.N0, self.N0 + self.Nitmax):
            self.Nit_loc = 0
            print("\n")
            print("Global iteration =  %s" % self.Nit)
            print("#" * 60)
            self.beta = 2 ** self.Nit
            # optimize it!

            try:
                popt = self.opt.optimize(p0)
            except nlopt.ForcedStop:
                print("forced stop...")
                # print(self.tot_obj_history)
                # print(self.Nit_loc)
                # print(self.tot_obj_history[-self.Nit_loc:])
                loc_obj = np.array(self.tot_obj_history[-self.Nit_loc :])
                loc_vars = np.array(self.param_history[-self.Nit_loc :])
                lb_ind = loc_vars >= self.opt.get_lower_bounds()
                ub_ind = loc_vars <= self.opt.get_upper_bounds()
                const_ind = np.all(lb_ind & ub_ind, axis=1)
                loc_obj = loc_obj[const_ind]
                loc_vars = loc_vars[const_ind]
                if self.typeopt is "max":
                    index = np.argmax(loc_obj)
                else:
                    index = np.argmin(loc_obj)

                # print(loc_obj[index])
                # print(loc_vars[index])
                popt = loc_vars[index]

            opt_f = self.opt.last_optimum_value()

            print(popt)
            # print("optimum at ", popt)
            print("-" * 45)
            print("   optimum   = ", opt_f)
            # print("result code = ", opt.last_optimize_result())
            result_code = self.opt.last_optimize_result()
            s0, s = self.opt_message(result_code)
            print(s0)
            print(s)
            # re-initialize for next global iteration
            p0 = np.copy(popt)
        return popt, opt_f, opt

    def opt_message(self, code):
        if code >= 0:
            s0 = "Successful termination"
        else:
            s0 = "Error"
        if code == 1:
            s = "Generic success return value."
        elif code == 2:
            s = "Optimization stopped because stopval was reached."
        elif code == 3:
            s = "Optimization stopped because ftol_rel or ftol_abs was reached."
        elif code == 4:
            s = "Optimization stopped because xtol_rel or xtol_abs was reached."
        elif code == 5:
            s = "Optimization stopped because maxeval was reached."
        elif code == 6:
            s = "Optimization stopped because maxtime was reached"
        elif code == -1:
            s = "Generic failure code."
        elif code == -2:
            s = "Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etc)."
        elif code == -3:
            s = "Ran out of memory."
        elif code == -4:
            s = "Halted because roundoff errors limited progress. (In this case, the optimization still typically returns a useful result.)"
        elif code == -5:
            s = "Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization’s nlopt_opt object opt from the user’s objective function or constraints."
        return s0, s

    def threshold_design(self, f_obj, p):
        print("\n")
        print("Final design")
        print("#" * 60)
        p_thres = np.zeros_like(p)
        threshold_val = self.p_interp
        threshold = np.linspace(0, 1, self.nthres + 1)
        for i in range(self.nthres):
            cond = np.logical_and(p >= threshold[i], p <= threshold[i + 1])
            p_thres[cond] = threshold_val[i]
        opt_thres = f_obj(p_thres, np.array([]), filt=False)
        print("-" * 45)
        print("   optimum   = ", opt_thres)
        return p_thres, opt_thres

    def get_threshold_design(self, p):
        p_thres = np.zeros_like(p)
        threshold_val = self.p_interp
        threshold = np.linspace(0, 1, self.nthres + 1)
        for i in range(self.nthres):
            cond = np.logical_and(p >= threshold[i], p <= threshold[i + 1])
            p_thres[cond] = threshold_val[i]
        return p_thres

    #
    # def grad_filter_param(self, p, xdes, ydes):
    #     if self.rfilt == 0:
    #         dpfilt = np.ones_like(p)
    #     else:
    #         dpfilt = []
    #         for xe, ye in zip(xdes, ydes):
    #             Ne = self.neighbourhood(xe, ye, xdes, ydes)
    #             xi, yi, pi = xdes[Ne], ydes[Ne], p[Ne]
    #             we = self.weight_filt(xe, ye, xe, ye)
    #             wi = self.weight_filt(xi, yi, xe, ye)
    #             dptmp = we / sum(wi)
    #             dpfilt.append(dptmp)
    #     # return np.array(dpfilt)
    #     return np.ones_like(p)

    # def grid2mesh(self, x_grid, y_grid, val_grid, xdes, ydes):
    #     f = sc.interpolate.interp2d(x_grid, y_grid, val_grid)
    #     i = 0
    #     x0 = np.zeros_like(xdes)
    #     for x, y in zip(xdes, ydes):
    #         x0[i] = f(x, y)
    #         i += 1
    #     return x0
    # def nodes2el(self, xnodes, ynodes, valnodes, xels, yels, method="linear"):
    #     points = np.vstack((xnodes, ynodes)).T
    #     valnodes.flags.writeable = True
    #     new_points = (xels, yels)
    #     re = sc.interpolate.griddata(points, valnodes.real, new_points, method=method)
    #     im = sc.interpolate.griddata(points, valnodes.imag, new_points, method=method)
    #     # re = self.grid2mesh(xnodes, ynodes, valnodes.real, xels, yels)
    #     # im = self.grid2mesh(xnodes, ynodes, valnodes.imag, xels, yels)
    #     return re + 1j * im
    # def proj(self, x):
    #     return x

    # def grad_simp(self, p):
    #     return self.m * (self.eps_max - self.eps_min) * p**(self.m - 1)
    #
    #
    # def grad_proj(self, x):
    #     return 0.5*self.beta/(np.tanh(self.thres*self.beta))*(1 - np.tanh(self.beta*(x - self.thres))**2)
