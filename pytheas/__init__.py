from .__about__ import (
    __version__,
    __author__,
    __author_email__,
    __copyright__,
    __website__,
    __license__,
    __status__,
    __description__,
)

__doc__ = __description__

from .basefem import BaseFEM
from .material import *
from .material.refractiveindex import Materials
from .material.genmat import MaterialDensity
from .periodic2D import *
from .scatt2D import *
from .scatt3D import *
from .periodic3D import *
from .banddiag import *

# from .optim import *
from .homogenization import *
from .tools import *
