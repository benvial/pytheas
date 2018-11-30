from .topopt import TopologyOptimization


def norm_vec(V):
    return np.sqrt(np.sum(V ** 2, axis=0))


class TopologyOptimization3D(TopologyOptimization):
    """A class for topology optimization in 3D

    """

    def __init__(self, fem):
        super().__init__(fem)

    n_x, n_y, n_z = 50, 50, 50

    @property
    def grid(self):
        xdes, ydes, zdes = self.fem.des[1].T
        x1 = np.linspace(min(xdes), max(xdes), self.n_x)
        y1 = np.linspace(min(ydes), max(ydes), self.n_y)
        z1 = np.linspace(min(zdes), max(zdes), self.n_z)
        return x1, y1, z1

    # Filter
    def neighbourhood(self, Xe, Xn):
        dX = np.array([Xe[0] - Xn[0, :], Xe[1] - Xn[1, :], Xe[2] - Xn[2, :]])
        return norm_vec(dX) <= self.rfilt

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

    def _filter_param(self, p):
        self.fem.print_progress("Filtering")
        Xdes = np.array(self.fem.des[1].T)
        if self.rfilt == 0:
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

    def make_xsym(self, qt, interp_method="cubic"):
        qt = self.mesh2grid(qt, interp_method=interp_method)
        qt = (qt + np.fliplr(qt)) / 2
        return self.grid2mesh(qt, interp_method=interp_method)

    def random_pattern(self, mat):
        im = mat.filtered_pattern[:, :, 0].T
        x_grid = np.linspace(self.grid[0][0], self.grid[0][-1], mat.n_x)
        y_grid = np.linspace(self.grid[1][0], self.grid[1][-1], mat.n_y)
        z_grid = np.linspace(self.grid[1][0], self.grid[1][-1], mat.n_z)
        x0 = self.grid2mesh(im, grid=(x_grid, y_grid, z_grid))
        return normalize(x0)

    def grid2mesh(self, val_grid, grid=None, interp_method="cubic"):
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

    def get_deq_deps(self):
        deq_deps = self.fem.get_deq_deps()
        return deq_deps

    def get_sensitivity(self, p, filt=True, proj=True):
        adjoint = self.get_adjoint()
        deq_deps = self.get_deq_deps()
        _, depsilon_dp = self.make_epsilon(p, filt=filt, proj=proj, grad=True)
        dotprod = np.sum(np.array([adjoint[i] * deq_deps[i] for i in range(3)]), axis=0)
        sens = self.dg_dp + 1 * np.real(dotprod * depsilon_dp)
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
