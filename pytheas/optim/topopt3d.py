from .topopt import TopOpt, normalize
import numpy as np
from scipy import interpolate


class TopOpt3D(TopOpt):
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

    def grid2mesh(self, val_grid, grid=None, interp_method="nearest"):
        des = self.fem.des[1].T
        if grid is None:
            grid = self.grid
        x0 = interpolate.interpn(grid, val_grid, des.T, method=interp_method)

        return np.array(x0).ravel()

    def mesh2grid(self, valdes, interp_method="nearest"):
        xdes, ydes, zdes = self.fem.des[1].T
        points = np.vstack((xdes, ydes, zdes)).T
        x_grid, y_grid, z_grid = self.grid
        valdes = np.require(valdes, requirements=["OWNDATA"])
        valdes.flags.writeable = True
        xg, yg, zg = np.meshgrid(x_grid, y_grid, z_grid)
        v = interpolate.griddata(
            points, valdes, (xg, yg, zg), method=interp_method, fill_value=0
        )
        v = np.transpose(v, [1, 0, 2])
        return v

    def get_deq_deps(self):
        return self.fem.get_deq_deps()

    def get_sensitivity(self, p, filt=True, proj=True):
        adjoint = self.get_adjoint()
        deq_deps = self.get_deq_deps()
        _, depsilon_dp = self.make_epsilon(p, filt=filt, proj=proj, grad=True)
        dotprod = np.sum(np.array([adjoint[i] * deq_deps[i] for i in range(3)]), axis=0)
        sens = self.dg_dp + 1 * np.real(dotprod * depsilon_dp)
        return sens

    def random_pattern(self, mat):
        # self=to
        im = mat.filtered_pattern
        # im = np.transpose(im, (1, 0, 2))
        x_grid = np.linspace(self.grid[0][0], self.grid[0][-1], mat.n_x)
        y_grid = np.linspace(self.grid[1][0], self.grid[1][-1], mat.n_y)
        z_grid = np.linspace(self.grid[2][0], self.grid[2][-1], mat.n_z)
        x0 = self.grid2mesh(im, grid=(x_grid, y_grid, z_grid))
        return normalize(x0)
