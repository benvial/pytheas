import warnings

warnings.filterwarnings("ignore", ".*GUI is implemented.*")
import numpy as np
import matplotlib

try:
    __IPYTHON__
except NameError:
    # print("Non interactive mode: setting matplotlib AGG backend")
    matplotlib.use("AGG")
else:
    try:
        import PyQt5
    except ImportError:
        # print("Interactive mode: setting matplotlib Qt4Agg backend")
        matplotlib.use("Qt4Agg")
    else:
        # print("Interactive mode: setting matplotlib Qt5Agg backend")
        matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator

plt.ion()

matplotlib.rcParams["xtick.major.size"] = 6
matplotlib.rcParams["xtick.major.width"] = 0.7
matplotlib.rcParams["xtick.minor.size"] = 0
matplotlib.rcParams["xtick.minor.width"] = 0
matplotlib.rcParams["ytick.major.size"] = 6
matplotlib.rcParams["ytick.major.width"] = 0.7
matplotlib.rcParams["ytick.minor.size"] = 0
matplotlib.rcParams["ytick.minor.width"] = 0
matplotlib.rcParams["savefig.dpi"] = 600
matplotlib.rcParams["lines.linewidth"] = 1.5
# matplotlib.rcParams['axes.linewidth'] = 0.7
matplotlib.rcParams["axes.labelsize"] = "medium"
matplotlib.rcParams["legend.fontsize"] = "small"
matplotlib.rcParams["font.size"] = 12
matplotlib.rcParams["axes.titlesize"] = "medium"
matplotlib.rcParams["xtick.labelsize"] = "x-small"
matplotlib.rcParams["ytick.labelsize"] = "x-small"

rc = {"axes.linewidth": 0.7}
sns.set_style("ticks", rc=rc)
# sns.set_context("poster", font_scale=1)

color1 = "#f4a688"
color2 = "#8d9db6"
color3 = "#27ae60"

aotomat_green = "#31a180"
aotomat_purple = "#a13152"
aotomat_yellow = "#ffe56e"
aotomat_orange = "#cc8800"
aotomat_light_gray = "#EEEEEE"
aotomat_gray_1 = "#BBBBBB"
aotomat_gray_2 = "#959393"
aotomat_gray_3 = "#888888"
aotomat_dark_gray_1 = "#333333"
aotomat_dark_gray_2 = "#444444"
aotomat_black = "#000000"


qmulLightGray = "#EEEEEE"
qmulGray1 = "#BBBBBB"
qmulGray2 = "#959393"
qmulGray3 = "#888888"
qmulDarkGray1 = "#333333"
qmulDarkGray2 = "#444444"
qmulBlack = "#000000"
qmulLightBlue = "#C1D3F3"
qmulBlue = "#1A428A"
qmulBlue1 = "#1D4998"
qmulBlue2 = "#1A4280"
qmulBlue3 = "#5A77AB"
qmulBlue4 = "#2157A5"
qmulRedPurple = "#C70151"
qmulLightorange = "#D36868"
qmulGreen = "#0D8175"

benchcolor1 = "#9E5B5B"  # %%% red %% frame title and rule
benchcolor2 = "#6FBD90"  # %%  green
benchcolor3 = "#424242"  # %%% dark grey/ black  %%normal text
benchcolor4 = "#03658C"  # %%   blue
benchcolor5 = "#E05A00"  # %% orange
benchcolor6 = "#424242"  # %%%% dark grey/ black
benchcolor7 = "#179D9D"  # %%%% green blue
benchcolor8 = "#B6855A"  # %%%% brown
benchcolor9 = "#9F61A8"  # %%%% purple
benchcolor10 = "#E09EDE"  # %%%% pink
benchcolor11 = "#A1A1A1"  # %%%% light gray
benchcolor12 = "#D0D0D0"  # %%%% lighter gray


cmap = sns.cubehelix_palette(8, start=0.5, rot=-0.75, as_cmap=True)


def colorbar(mappable, **kwargs):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax, **kwargs)


def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1.0, N), (0.0, 0.0, 0.0, 0.0)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1.0, N + 1)
    cdict = {}
    for ki, key in enumerate(("red", "green", "blue")):
        cdict[key] = [
            (indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
            for i in range(N + 1)
        ]
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, 1024)
