import numpy as np
import nlopt
from pytheas.tools.plottools import *
import numpy as np
import scipy as sc

# from scipy.interpolate import splev, splrep
from scipy.interpolate import PchipInterpolator  # mono_cubic_interp


class TopologyOptimization:
    """A class for topology optimization

    """

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
    stopval = 1.0
    maxeval = 50  # maximum of function evaluation
    Nitmax = 8  # maximum number of global iterations
    Nit = 0  # initialize global iteration number
    beta = 2 ** Nit  # projection parameter
    rfilt = 0  # filter radius
    plotconv = False
    dg_dp = 0
    eps_min, eps_max = 1, 5
    nthres = 2
    beta = 1
    log_opt = False
    dp = 1e-3

    @property
    def thres(self):
        return np.linspace(0, 1, self.nthres + 1)[1 : (self.nthres)]

    @property
    def p_interp(self):
        return np.linspace(0, 1, self.nthres)

    def simp(self, p):
        # k = min(self.nthres - 1, 3)
        if self.nthres == 2:
            out = (
                self.eps_interp[1] - self.eps_interp[0]
            ) * p ** self.m + self.eps_interp[0]
        else:
            # tsimp_re = splrep(self.p_interp, self.eps_interp.real, k=k)
            # tsimp_im = splrep(self.p_interp, self.eps_interp.imag, k=k)
            # out = splev(p, tsimp_re) + 1j*splev(p, tsimp_im)
            tsimp_re = PchipInterpolator(self.p_interp, self.eps_interp.real)
            tsimp_im = PchipInterpolator(self.p_interp, self.eps_interp.imag)
            out = tsimp_re(p ** self.m) + 1j * tsimp_im(p ** self.m)
        return out

    def proj(self, x):
        p = 0
        for thres in self.thres:
            ptmp = 0.5 + 0.5 * np.tanh(self.beta * (x - thres)) / (
                np.tanh(thres * self.beta)
            )
            p += self.normalize(ptmp)
            # p += ptmp #self.normalize(ptmp)
        # return p / (self.nthres - 1)
        return x

    # def proj(self, x):
    #     return x

    def grad_simp(self, p):
        f = self.simp(p)
        df = self.simp(p + self.dp)
        return (df - f) / self.dp

    def grad_proj(self, p):
        f = self.proj(p)
        df = self.proj(p + self.dp)
        return (df - f) / self.dp

    # def grad_simp(self, p):
    #     return self.m * (self.eps_max - self.eps_min) * p**(self.m - 1)
    #
    #
    # def grad_proj(self, x):
    #     return 0.5*self.beta/(np.tanh(self.thres*self.beta))*(1 - np.tanh(self.beta*(x - self.thres))**2)

    # Filter
    def neighbourhood(self, xe, ye, xn, yn):
        D = np.sqrt((xn - xe) ** 2 + (yn - ye) ** 2)
        return D <= self.rfilt

    def weight_filt(self, x, y, xe, ye):
        D = np.sqrt((x - xe) ** 2 + (y - ye) ** 2)
        s = self.rfilt / 2
        return np.exp(-0.5 * (D / s) ** 2)

    def filter_param(self, p, xdes, ydes):
        if self.rfilt == 0:
            pfilt = p
        else:
            pfilt = []
            for xe, ye, pe in zip(xdes, ydes, p):
                Ne = self.neighbourhood(xe, ye, xdes, ydes)
                xi, yi, pi = xdes[Ne], ydes[Ne], p[Ne]
                wi = self.weight_filt(xi, yi, xe, ye)
                ptmp = sum(wi * pi) / sum(wi)
                pfilt.append(ptmp)
        return np.array(pfilt)
        # return np.array(p)

    def grad_filter_param(self, p, xdes, ydes):
        if self.rfilt == 0:
            dpfilt = np.ones_like(p)
        else:
            dpfilt = []
            for xe, ye, pe in zip(xdes, ydes, p):
                Ne = self.neighbourhood(xe, ye, xdes, ydes)
                xi, yi, pi = xdes[Ne], ydes[Ne], p[Ne]
                we = self.weight_filt(xe, ye, xe, ye)
                wi = self.weight_filt(xi, yi, xe, ye)
                dptmp = we / sum(wi)
                dpfilt.append(dptmp)
        return np.array(dpfilt)
        # return np.ones_like(p)

    def make_design_rect_grid(self, xdes, ydes, nx, ny):
        x1 = np.linspace(min(xdes), max(xdes), nx)
        y1 = np.linspace(min(ydes), max(ydes), ny)
        return x1, y1

    def normalize(self, x0):
        if x0.min() == x0.max():
            return x0
        else:
            return (x0 - x0.min()) / (x0.max() - x0.min())

    def random_pattern(self, xdes, ydes, x_grid, y_grid, mat):
        im = mat.filtered_pattern[:, :, 0].T
        x0 = self.grid2mesh(x_grid, y_grid, im, xdes, ydes)
        return self.normalize(x0)

    def grid2mesh(self, x_grid, y_grid, val_grid, xdes, ydes):
        f = sc.interpolate.interp2d(x_grid, y_grid, val_grid)
        x0 = [f(x, y) for x, y in zip(xdes, ydes)]
        return np.array(x0).ravel()

    #
    # def grid2mesh(self, x_grid, y_grid, val_grid, xdes, ydes):
    #     f = sc.interpolate.interp2d(x_grid, y_grid, val_grid)
    #     i = 0
    #     x0 = np.zeros_like(xdes)
    #     for x, y in zip(xdes, ydes):
    #         x0[i] = f(x, y)
    #         i += 1
    #     return x0

    def mesh2grid(self, xdes, ydes, valdes, x_grid, y_grid, method="nearest"):
        points = np.vstack((xdes, ydes)).T
        valdes.flags.writeable = True
        xg, yg = np.meshgrid(x_grid, y_grid)
        v = sc.interpolate.griddata(points, valdes, (xg, yg), method=method)
        return v

    # def nodes2el(self, xnodes, ynodes, valnodes, xels, yels, method="linear"):
    #     points = np.vstack((xnodes, ynodes)).T
    #     valnodes.flags.writeable = True
    #     new_points = (xels, yels)
    #     re = sc.interpolate.griddata(points, valnodes.real, new_points, method=method)
    #     im = sc.interpolate.griddata(points, valnodes.imag, new_points, method=method)
    #     # re = self.grid2mesh(xnodes, ynodes, valnodes.real, xels, yels)
    #     # im = self.grid2mesh(xnodes, ynodes, valnodes.imag, xels, yels)
    #     return re + 1j * im

    def make_epsilon(self, p, xdes, ydes, filt=True):
        p_filt = self.filter_param(p, xdes, ydes)
        p_proj = self.proj(p_filt)
        epsilon = self.simp(p_proj)
        return epsilon

    def depsilon_dp(self, p, xdes, ydes):
        p_filt = self.filter_param(p, xdes, ydes)
        p_proj = self.proj(p_filt)
        epsilon = self.simp(p_proj)
        dpfilt_dp = self.grad_filter_param(p, xdes, ydes)
        dpproj_dpfilt = self.grad_proj(p)
        depsilon_dpproj = self.grad_simp(p)
        return depsilon_dpproj * dpproj_dpfilt * dpfilt_dp

    def get_objective(self, fem):
        g = fem.get_objective()
        if self.log_opt:
            goal = np.log(g)
        else:
            goal = g
        return goal

    def sensitivity(self, p, xdes, ydes, adjoint, deq_deps):
        sens = self.dg_dp + adjoint * deq_deps * self.depsilon_dp(p, xdes, ydes)
        if self.log_opt:
            p_min = 1e-3
            p_thres = np.copy(p)
            p_thres[p <= p_min] = p_min
            sens /= p_thres
        return sens

    def plot_design(
        self, ax, xdes, ydes, varplot, x_grid, y_grid, typeplot="interp", **kwargs
    ):
        if typeplot is "tri":
            triang = matplotlib.tri.Triangulation(xdes, ydes)
            xmid = xdes[triang.triangles].mean(axis=1)
            ymid = ydes[triang.triangles].mean(axis=1)
            cf = ax.tripcolor(triang, varplot, shading="flat")
        elif typeplot is "interp":
            varplot = self.mesh2grid(xdes, ydes, varplot, x_grid, y_grid, **kwargs)
            cf = ax.imshow(varplot)
        else:
            raise TypeError("Wrong typeplot specified: choose between interp and tri")
        cbar = plt.colorbar(cf, fraction=0.046, pad=0.04)
        # ax.set_title(r'Permittivity $\varepsilon$ in the design domain')
        ax.axis("image")
        ax.axis("off")

    def plot_convergence(self, ax, OBJ):
        if len(OBJ) > 20:
            styleplot = "-"
        else:
            styleplot = "-o"
        ax.plot(OBJ, styleplot, color=aotomat_green)
        # x_obj = list(range(len(OBJ)))
        # ax.fill_between(x_obj, OBJ, 0., color=aotomat_green, alpha=0.2)
        ax.set_title(
            "Convergence: "
            + "global iteration = "
            + str(self.Nit)
            + ", current iteration = "
            + str(len(OBJ) - 1)
        )
        ax.set_xlabel(r"iteration number $N$")
        ax.set_ylabel(r"objective $\varphi$")
        # ax4.grid(True)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.axis("tight")

    def plot_while_solving(
        self, varplot, xdes, ydes, x_grid, y_grid, OBJ, title="", **kwargs
    ):
        # print("beta = ", beta)
        plt.clf()
        ax1 = plt.subplot(211, aspect="equal")
        self.plot_design(ax1, xdes, ydes, varplot, x_grid, y_grid, **kwargs)
        ax1.set_title(title)
        ax2 = plt.subplot(212)
        self.plot_convergence(ax2, OBJ)
        ax2.set_xlim((0, len(OBJ) - 1))
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
        opt.set_stopval(self.stopval)
        opt.set_xtol_rel(self.ptol_rel)
        opt.set_ftol_rel(self.ftol_rel)
        ################################################################
        ############### OPTIMIZATION - global iterations ###############
        ################################################################
        for self.Nit in range(self.Nitmax):
            print("\n")
            print("Global iteration =  %s" % self.Nit)
            print("#" * 60)
            self.beta = 2 ** self.Nit

            # optimize it!
            popt = opt.optimize(p0)
            opt_f = opt.last_optimum_value()
            # print("optimum at ", popt)
            print("-" * 45)
            print("   optimum   = ", opt_f)
            # print("result code = ", opt.last_optimize_result())
            result_code = opt.last_optimize_result()
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
