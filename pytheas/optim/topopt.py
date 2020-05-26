import numpy as np
import nlopt
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from matplotlib.tri import Triangulation
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
from ..material import genmat

# from pytheas.tools.plottools import *


def opt_message(code):
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


def normalize(x0):
    if x0.min() == x0.max():
        return x0
    else:
        return (x0 - x0.min()) / (x0.max() - x0.min())


def norm_vec(V):
    return np.sqrt(np.sum(V ** 2, axis=0))


def simp_(p, val, m=1):
    nthres = len(val)
    if nthres == 2:
        out = (val[1] - val[0]) * p ** m + val[0]
    else:
        p_interp = np.linspace(0, 1, nthres)
        tsimp_re = interpolate.PchipInterpolator(p_interp, val.real)
        tsimp_im = interpolate.PchipInterpolator(p_interp, val.imag)
        out = tsimp_re(p ** m) + 1j * tsimp_im(p ** m)
    return out


def neighbourhood_(Xe, Xn, rfilt):
    dX_ = []
    for xe, xn in zip(Xe, Xn):
        dX_.append(xe - xn)
    dX = np.array(dX_)
    return norm_vec(dX) <= rfilt


class TopOpt:
    """A class for topology optimization"""

    def __init__(self, fem):
        self.obj_history = []
        self.param_history = []
        self.tot_obj_history = []
        self.fem = fem
        self.verbose = 1

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
    dp = 1e-11
    n_x, n_y, n_z = 100, 100, 1
    force_xsym = False
    filter_param_function = "new"

    callback = None

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
        xdes, ydes, _ = self.fem.des[1].T
        x1 = np.linspace(min(xdes), max(xdes), self.n_x)
        y1 = np.linspace(min(ydes), max(ydes), self.n_y)
        return x1, y1

    def _simp(self, p):
        return simp_(p, self.eps_interp, m=self.m)

    def simp(self, p, grad=True):
        g = None
        f = self._simp(p)
        if grad:
            df = self._simp(p + self.dp)
            g = (df - f) / self.dp
        return f, g

    def _proj(self, x):
        p = 0
        for thres in self.thres:
            ptmp = np.tanh(self.beta * (x - thres)) / (np.tanh(thres * self.beta))
            p += ptmp
        p = (p + self.nthres - 1) / (2 * self.nthres - 2)
        return p

    def proj(self, p, grad=True):
        g = None
        f = self._proj(p)
        if grad:
            df = self._proj(p + self.dp)
            g = (df - f) / self.dp
        return f, g

    # Filter
    def neighbourhood(self, Xe, Xn):
        return neighbourhood_(Xe, Xn, self.rfilt)

    def weight_filt(self, Xe, Xn):
        dX = np.array([Xe[0] - Xn[0, :], Xe[1] - Xn[1, :], Xe[2] - Xn[2, :]])
        D = norm_vec(dX)
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

    def old_filter_param(self, p):
        self.fem._print_progress("Filtering")
        Xdes = np.array(self.fem.des[1].T)
        isfilt = np.array(self.rfilt) == 0
        nofilt = isfilt.all()
        if nofilt:
            pfilt = p
        else:
            pfilt = np.zeros_like(p)
            for i, Xe in enumerate(Xdes.T):
                Ne = self.neighbourhood(Xe, Xdes)
                Xneig = Xdes[:, Ne]
                pneig = p[Ne]
                w = self.weight_filt(Xe, Xneig)
                pfilt[i] = np.sum(w * pneig) / np.sum(w)
        return pfilt

    def new_filter_param(self, p):
        self.fem._print_progress("Filtering")
        isfilt = np.array(self.rfilt) == 0
        nofilt = isfilt.all()
        if nofilt:
            pfilt = p
        else:
            if hasattr(self.rfilt, "__len__"):
                if len(self.rfilt) == 2:
                    sigma = np.array(
                        [
                            np.sqrt(self.rfilt[0] * self.n_x),
                            np.sqrt(self.rfilt[1] * self.n_y),
                        ]
                    )
                else:
                    if len(self.rfilt) == 3:
                        sigma = np.array(
                            [
                                np.sqrt(self.rfilt[0] * self.n_x),
                                np.sqrt(self.rfilt[1] * self.n_y),
                                np.sqrt(self.rfilt[2] * self.n_z),
                            ]
                        )
            else:
                if self.fem.dim == 2:
                    sigma = np.sqrt(self.rfilt * np.array([self.n_x, self.n_y]))
                else:
                    if self.fem.dim == 3:
                        sigma = np.sqrt(
                            self.rfilt * np.array([self.n_x, self.n_y, self.n_z])
                        )

            p_grid = self.mesh2grid(p)
            p_filt_scipy = gaussian_filter(p_grid, sigma=np.flipud(sigma))
            pfilt = self.grid2mesh(p_filt_scipy)
        return pfilt

    def _filter_param(self, p):
        if self.filter_param_function == "old":
            return self.old_filter_param(p)
        else:
            if self.filter_param_function == "new":
                return self.new_filter_param(p)
            else:
                raise TypeError(
                    "Wrong filter_param_function specified: choose between new and old"
                )

    def filter_param(self, p, grad=True):
        g = None
        f = self._filter_param(p)
        if grad:
            df = self._filter_param(p + self.dp)
            g = (df - f) / self.dp
        return f, g

    def make_xsym(self, qt, interp_method="cubic"):
        qt = self.mesh2grid(qt, interp_method="nearest")
        qt = (qt + np.fliplr(qt)) / 2
        return self.grid2mesh(qt, interp_method=interp_method)

    def random_pattern(self, mat, discrete=False):
        im = mat.filtered_pattern[:, :, 0].T
        x_grid = np.linspace(self.grid[0][0], self.grid[0][-1], mat.n_x)
        y_grid = np.linspace(self.grid[1][0], self.grid[1][-1], mat.n_y)
        im_mesh = self.grid2mesh(im, grid=(x_grid, y_grid))
        im_mesh = normalize(im_mesh)
        if discrete:
            return genmat.make_discrete_pattern(im_mesh, mat.threshold_val)
        else:
            return im_mesh

    def grid2mesh(self, val_grid, grid=None, interp_method="cubic"):
        xdes, ydes, zdes = self.fem.des[1].T
        if grid:
            x_grid, y_grid = grid
        else:
            x_grid, y_grid = self.grid
        f = interpolate.interp2d(x_grid, y_grid, val_grid, kind=interp_method)
        x0 = [f(x, y) for x, y in zip(xdes, ydes)]
        return np.array(x0).ravel()

    def mesh2grid(self, valdes, interp_method="cubic"):
        xdes, ydes, zdes = self.fem.des[1].T
        points = np.vstack((xdes, ydes)).T
        x_grid, y_grid = self.grid
        valdes = np.require(valdes, requirements=["OWNDATA"])
        valdes.flags.writeable = True
        xg, yg = np.meshgrid(x_grid, y_grid)
        v = interpolate.griddata(
            points, valdes, (xg, yg), method=interp_method, fill_value=0
        )
        return v

    def make_epsilon(self, p, filt=True, proj=True, grad=False):
        self.fem._print_progress("Building permittivity")
        p_filt, dp_filt_dp = (
            self.filter_param(p, grad=grad) if filt else (p, np.ones_like(p))
        )
        p_proj, dp_proj_dpfilt = (
            self.proj(p_filt, grad=grad) if proj else (p_filt, np.ones_like(p))
        )
        epsilon, depsilon_dpproj = self.simp(p_proj, grad=grad)
        if grad:
            return epsilon, depsilon_dpproj * dp_proj_dpfilt * dp_filt_dp
        else:
            return epsilon

    def get_objective(self):
        goal = self.fem.get_objective()
        return goal

    def get_adjoint(self):
        return self.fem.get_adjoint()

    def get_deq_deps(self, interp_method="nearest"):
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
        self.fem._print_progress("Retrieving sensitivity")
        adjoint = self.get_adjoint()
        deq_deps = self.get_deq_deps(interp_method=interp_method)
        _, depsilon_dp = self.make_epsilon(p, filt=filt, proj=proj, grad=True)
        sens = self.dg_dp + np.real(adjoint * deq_deps * depsilon_dp)
        return sens

    def plot_design(
        self, ax, varplot, typeplot="interp", cmap=None, extent=None, **kwargs
    ):
        if typeplot is "tri":
            xdes, ydes, _ = self.fem.des[1].T
            triang = Triangulation(xdes, ydes)
            cf = ax.tripcolor(triang, varplot, shading="flat", cmap=cmap)
        elif typeplot is "interp":
            varplot = self.mesh2grid(varplot, **kwargs)
            cf = ax.imshow(np.flipud(varplot), cmap=cmap, extent=extent)
        else:
            raise TypeError("Wrong typeplot specified: choose between interp and tri")
        plt.colorbar(cf, fraction=0.046, pad=0.04)
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
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.axis("tight")

    def plot_while_solving(self, varplot, title="", **kwargs):
        plt.clf()
        ax1 = plt.subplot(211, aspect="equal")
        self.plot_design(ax1, varplot, **kwargs)
        ax1.set_title(title)
        ax2 = plt.subplot(212)
        self.plot_convergence(ax2)
        if self.Nit_tot > 1:
            ax2.set_xlim((0, self.Nit_tot - 1))
        plt.tight_layout()
        plt.show()
        plt.pause(0.01)

    def main_loop_topopt(self, f_obj, p0):
        nvar = len(p0)
        if self.verbose:
            print("\n")
            print("#" * 60)
            print("Topology optimization with %s variables" % nvar)
            print("#" * 60)
        # instanciation opt object, optimization algorithm
        opt = nlopt.opt(self.algorithm, nvar)

        ################################################################
        ############### OPTIMIZATION - global iterations ###############
        ################################################################
        for self.Nit in range(self.N0, self.N0 + self.Nitmax):

            self.beta = 2 ** self.Nit
            opt.set_lower_bounds(self.pmin)
            opt.set_upper_bounds(self.pmax)
            if self.typeopt is "max":
                opt.set_max_objective(f_obj)
            elif self.typeopt is "min":
                opt.set_min_objective(f_obj)
            else:
                raise TypeError("Wrong typeopt specified: choose between max and min")
            opt.set_maxeval(self.maxeval)
            if self.stopval:
                opt.set_stopval(self.stopval)
            opt.set_xtol_rel(self.ptol_rel)
            opt.set_ftol_rel(self.ftol_rel)
            self.opt = opt
            self.Nit_loc = 0
            if self.verbose:
                print("\n")
                print("Global iteration =  %s" % self.Nit)
                print("#" * 60)
            # optimize it!
            try:
                popt = self.opt.optimize(p0)
            except nlopt.ForcedStop:
                if self.verbose:
                    print("Forced stop...")
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
                popt = loc_vars[index]

            opt_f = self.opt.last_optimum_value()
            result_code = self.opt.last_optimize_result()
            s0, s = opt_message(result_code)
            if self.verbose:
                print(popt)
                print("-" * 45)
                print("   optimum   = ", opt_f)
                print(s0)
                print(s)
            # re-initialize for next global iteration
            p0 = np.copy(popt)
            if self.callback:
                self.callback(self)
        return popt, opt_f, opt

    def threshold_design(self, f_obj, p):
        p_thres = self.get_threshold_design(p)
        opt_thres = f_obj(p_thres, np.array([]), filt=False)
        if self.verbose:
            print("\n")
            print("Final design")
            print("#" * 60)
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
