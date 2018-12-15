from .topopt import TopOpt
import numpy as np


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

    def get_deq_deps(self):
        return self.fem.get_deq_deps()

    def get_sensitivity(self, p, filt=True, proj=True):
        adjoint = self.get_adjoint()
        deq_deps = self.get_deq_deps()
        _, depsilon_dp = self.make_epsilon(p, filt=filt, proj=proj, grad=True)
        dotprod = np.sum(np.array([adjoint[i] * deq_deps[i] for i in range(3)]), axis=0)
        sens = self.dg_dp + 1 * np.real(dotprod * depsilon_dp)
        return sens
