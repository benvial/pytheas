import numpy as np
from scipy import ndimage
pi = np.pi

class MaterialDensity:
    """A class for generating material densities"""

    n_x = 2**8
    n_y = 2**8
    n_z = 2**8
    nb_threshold = 2
    ratio_filter = 8
    xsym = False
    ysym = False
    zsym = False
    indep_z = False
    sym8 = False
    _threshold_val = np.linspace(0, 1, nb_threshold)

    @property
    def mat_grid(self):
        X = np.arange(0, self.n_x, 1)
        Y = np.arange(0, self.n_y, 1)
        Z = np.arange(0, self.n_z, 1)
        return np.meshgrid(X, Y, Z, indexing='ij')

    @property
    def mat_rand(self):
        return np.random.random((self.n_x, self.n_y, self.n_z))

    @property
    def pattern_rand(self):
        p = self.p_seed
        return self.make_sym(p)

    def make_sym(self, p):
        if not self.sym8:
            if self.xsym:
                p[-int(self.n_x / 2):, :,
                  :] = np.flipud(p[0:int(self.n_x / 2), :, :])
            if self.ysym:
                p[:, -int(self.n_y / 2):,
                  :] = np.fliplr(p[:, 0:int(self.n_y / 2), :])

        if self.zsym:
            p[:, -int(self.n_y / 2):,
              :] = np.fliplr(p[:, 0:int(self.n_y / 2), :])
        if self.sym8:
            assert self.n_x == self.n_y, "x and y sampling must be the same!"
            p1 = p[0:int(self.n_x / 2), 0:int(self.n_y / 2), :]
            p1 = 0.5 * (p1 + np.transpose(p1, axes=[1, 0, 2]))
            p[0:int(self.n_x / 2), 0:int(self.n_y / 2), :] = p1
            p[0:int(self.n_x / 2), -int(self.n_y / 2):, :] = np.fliplr(p1)
            p[-int(self.n_x / 2):, 0:int(self.n_y / 2), :] = np.flipud(p1)
            p[-int(self.n_x / 2):, -int(self.n_y / 2)
                   :, :] = np.flipud(np.fliplr(p1))
        if self.indep_z:
            p = np.dstack([p[:, :, 0]] * self.n_z)

        return p

    @property
    def sigma(self):
        return np.array([self.n_x, self.n_y, self.n_z]) / np.array(self.ratio_filter)

    # @property
    # def pattern_cont(self):
    #     return self.make_continuous_pattern()
    #
    # @property
    # def pattern_disc(self):
    #     return self.make_discrete_pattern(self.pattern_cont)

    @property
    def filtered_pattern(self):
        # apply gaussian filter
        return filter_pattern(self.pattern_rand, self.sigma)

    @property
    def normalized_pattern(self):
        # normalize between 0 and 1
        return normalize(self.filtered_pattern)

    @property
    def threshold(self):
        return np.linspace(0, 1, self.nb_threshold + 1)

    # @property
    # def _threshold_val(self):
    #     return np.linspace(0, 1, self.nb_threshold)

    @property
    def threshold_val(self):
        return np.linspace(0, 1, self.nb_threshold)

    @threshold_val.setter
    def threshold_val(self, value):
        self._threshold_val = value

    @property
    def discrete_pattern(self):
        return make_discrete_pattern(self.normalized_pattern, self._threshold_val)
    #

    def plot_pattern(self, fig, ax, extent=None, interpolation=None, cmap="Blues", indexz=0):
        mask = self.discrete_pattern[:, :, indexz]
        im = ax.imshow(np.fliplr(mask).T, interpolation=interpolation,
                       cmap=cmap, vmin=0, vmax=1, extent=extent)
        ax.set_axis_off()
        ax.set_title('material distribution')
        fig.colorbar(im, fraction=0.046, pad=0.04)


    def fourrier_pattern(self, c, scale=(1,1), norma=True):
        Xg, Yg, Zg = self.mat_grid
        Xn = 1 * Xg / (self.n_x - 1) - 1 / 2
        Yn = 1 * Yg / (self.n_y - 1) - 1 / 2
        PAT = 0
        Nhx, Nhy = c.shape
        for ix in range(Nhx):
            for iy in range(Nhy):
                kx = 2 * pi * (-Nhx/2 + 1/2 + ix)
                ky = 2 * pi * (-Nhy/2 + 1/2 + iy)
                xsi = scale[0] * kx * Xn + scale[1] * ky * Yn
                p = c[ix, iy] * np.exp(1j*xsi)
                PAT += p
        # # PAT = PAT[:, :, 0]
        # if self.sym8:
        #     assert self.n_x == self.n_y, "x and y sampling must be the same!"
        #     PAT = 0.5*(PAT  + (np.fliplr(PAT)))
        #     PAT = 0.5*(PAT  + (np.flipud(PAT)))
        #     PAT = 0.5*(PAT  + (PAT).T)
        # elif self.xsym:
        #     PAT = 0.5*(PAT  + (np.fliplr(PAT)))
        #     if self.ysym:
        #         PAT = 0.5*(PAT  + (np.flipud(PAT)))
        # elif self.ysym:
        #     PAT = 0.5*(PAT  + (np.flipud(PAT)))
        # else:
        #     pass
        # PAT = PAT/np.sqrt((Nhx-1)/2*(Nhy-1)/2)
        PAT = np.real(PAT)
        PAT = self.make_sym(PAT)
        if norma:
            PAT = normalize(PAT)
        PAT = make_discrete_pattern(PAT, self._threshold_val)
        return PAT





def multi_period(pat, npx, npy):
    patper = pat
    for ix in range(npx-1):
        patper = np.vstack((patper, pat))
    patper1 = patper
    for iy in range(npy-1):
        patper1 = np.hstack((patper1, patper))

    return patper1

def normalize(im):
    return (im - im.min()) / (im.max() - im.min())


def filter_pattern(p, sigma):
    return ndimage.gaussian_filter(p, sigma, order=0)



def make_discrete_pattern(im, threshold_val):
    mask = np.zeros(im.shape)
    L = len(threshold_val)
    threshold = np.linspace(0, 1, L + 1)
    for i in range(L):
        cond = np.logical_and(im > threshold[i], im <= threshold[i + 1])
        mask[cond] = threshold_val[i]
    return mask


def rot_z(t):
    c = np.cos(t)
    s = np.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def ell_shapes(mat_grid, rloc=[0, 0, 0], rwidth=[0.1, 0.1, 0.1], m=2):
    coords1 = 0
    N = np.shape(mat_grid)[1:4]
    for i in range(3):
        rcoords = 0
        if N[i] is not 1:
            rcoords = mat_grid[i] / (N[i] - 1) - 0.5
        coords1 += ((rcoords - rloc[i]) / (1 * rwidth[i]))**m
    # return np.exp(-sum(coords1))
    # return np.exp(-coords1)
    return coords1 < 1


def ell_shapes_array(Nb, mat_grid, rloc, rwidth, m=2):
    b = 0
    for n in range(Nb):
        rloc1 = rloc[n, :]
        rwidth1 = rwidth[n, :]
        btmp = ell_shapes(mat_grid, rloc=rloc1, rwidth=rwidth1, m=m)
        b += btmp
        # if (b == 2).any():
        #     b -= btmp
    b[b != 0] = 1
    return b


def random_fibers(f, mat_grid, d0=0.1, dr_ratio=0.3,
                  touching=True, touching_ratio=0.01, m=2):
    b = 0
    n = 0
    F = -1
    # while np.abs(F - f)>1e-2:
    i0 = 0
    i1 = 0
    btmp_list = []
    tol_f = 5e-3
    cond_stop = True
    while cond_stop:
        rloc = np.random.rand(3) - 0.5
        rloc[2] = 0
        r0 = d0 / 2
        rwidth = r0 + 2 * (0.5 - 1 * np.random.rand(3)) * dr_ratio * r0
        rwidth[1] = rwidth[0]
        rloc = rloc * (1 - 2 * max(rwidth)) * 0.99
        btmp = ell_shapes(mat_grid, rloc=rloc, rwidth=rwidth, m=m)
        btmp_list.append(btmp)
        b += btmp
        if touching:
            bt = np.copy(b)
            bt[b == 1] = 0
            bt[b == 2] = 1
            t = np.mean(np.mean(bt))
            # print("touching_ratio = ", t)

            if t > touching_ratio * r0:
                # b -= np.sum(btmp_list[-i1:], axis=0)
                b -= btmp
                i1 += 1

        else:
            if (b == 2).any():
                b -= btmp

        bt = np.copy(b)
        bt[bt != 0] = 1
        F = np.mean(np.mean(bt))
        d_f = np.abs(F - f)

        if F >= (f + tol_f):
            b -= btmp
            # b -= np.sum(btmp_list[-i1:], axis=0)
            i1 += 1
        i0 += 1
        # print("i0 = ", i0)
        # print("i1 = ", i1)
        if i0 > 10000:
            print("restarting...")
            print(F)
            b = np.zeros_like(b)
            i0 = 0
            i1 = 0

        d_f = np.abs(F - f)
        # print(F)
        # print(d_f)
        cond_stop = (d_f > tol_f)
        # print(cond_stop)
    b[b != 0] = 1
    return b, F


if __name__ == '__main__':
    from pytheas.tools.plottools import *
    # plt.close("all")
    plt.clf()
    np.random.seed(12)

    mat = MaterialDensity()
    mat.n_x, mat.n_y, mat.n_z = 2**8, 2**8, 1

    mat.p_seed = np.random.random((mat.n_x, mat.n_y, mat.n_z))
    mat.nb_threshold = 4
    # mat.xsym = True
    # mat.ysym = True
    mat.sym8 = True
    mat.ratio_filter = [11, 11, 1]
    #
    # Xg, Yg, Zg = mat.mat_grid
    # sx, sy = mat.n_x / 2.2, mat.n_y / 2.2
    #
    # N = 8
    # Xa = (Xg - mat.n_x / 2)**N / sx**N
    # Ya = (Yg - mat.n_y / 2)**N / sy**N
    # U = 1
    # mat.p_seed = U
    #
    # plt.figure()
    # rloc = [0, 0, 0]
    # rwidth = [0.1, 0.2, 0.3]
    # m = 2
    # b = ell_shapes(mat.mat_grid, rloc=rloc, rwidth=rwidth, m=m)
    # # b = ndimage.interpolation.rotate(b, 30, axes=(1, 2), order=0, reshape=False)
    # # b = b[:, :, int(mat.n_z/2)]
    # # b = b[int(mat.n_x/2), :, :]
    # b = b[:, :, 0]
    # # b = normalize(b)
    # # b = make_discrete_pattern(b, mat.nb_threshold)
    # #
    # plt.imshow(b)
    # plt.colorbar()
    #
    # f = 0.3
    # b = random_fibers(f, mat.mat_grid, touching=True)[:, :, 0]
    #
    # # b = normalize(b)
    # # b = make_discrete_pattern(b, mat.nb_threshold)
    # plt.clf()
    # plt.imshow(b)
    # plt.colorbar()
    # csd
    #
    # # genmat.make_discrete_pattern(discrete_pattern)
    #
    # #
    # # plt.imshow( U[:,:,0])
    # # plt.colorbar()
    # # cd
    #
    # F = np.exp(-Xa) * np.exp(-Ya)
    # # F = np.zeros_like(mat.p_seed)
    # ncx = int(mat.n_x / 17)
    # # F[ncx:-ncx, ncx:-ncx,:] = 1
    # U = mat.normalized_pattern * F
    # # U = mat.mat_rand*F
    # # Udisc = filter_pattern(Udisc)
    # # Udisc = genmat.make_discrete_pattern(U, mat.nb_threshold)
    # sigma = mat.sigma * 0.01
    # Udisc = filter_pattern(U, sigma)
    # Udisc = normalize(Udisc)
    # Udisc = make_discrete_pattern(Udisc, mat.nb_threshold)
    #
    # fig, ax = plt.subplots(2)
    # ax[0].imshow(U[:, :, 0])
    #
    # ax[1].imshow(Udisc[:, :, 0])
    #
    # fig, ax = plt.subplots(1)
    # mat.plot_pattern(fig, ax, cmap=cmap)
    #
    # xsa

    np.random.seed(3)

    # mat.threshold_val = (mat.threshold_val)
    # mat.threshold_val = np.random.permutation(mat.threshold_val)
    iplot = range(33)
    for i in iplot:
        # mat.threshold_val = np.random.permutation(mat.threshold_val)
        plt.clf()
        fig = plt.gcf()
        ax = plt.gca()
        mat.plot_pattern(fig, ax, cmap=cmap)
        plt.draw()
        plt.pause(0.1)
