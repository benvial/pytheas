# -*- coding: utf-8 -*-
"""
Importing refractive index from a database
==========================================

Retrieve and plot the refractive index of a material in the refractiveindex.info 
data.
"""

# Code source: Benjamin Vial
# License: MIT

from pytheas.material.refractiveindex import *
from pytheas.tools.plottools import *

##############################################################################
# We can get the refractive index from tabulated data or a formula using the
# database in the :py:mod:`pytheas.material` module.
# We will import the measured data from the reference `Johnson and Christy`_ [JC1972]_.
# We first specify the file :code:`yamlFile` we want to import:

yamlFile = "main/Au/Johnson.yml"


##############################################################################
# We then get the wavelength bounds from the data (in microns) and create a
# wavelength range to interpolate:

bounds = getRange(yamlFile)
lambdas = np.linspace(bounds[0], bounds[1], 300)

##############################################################################
# Then get the refractive index data:

ncomplex = get_complex_index(lambdas, yamlFile)
epsilon = (ncomplex**2)

##############################################################################
# And finally plot it:

plt.close('all')
fig, ax = plt.subplots(1, figsize=(6,4))
plt.plot(lambdas , epsilon.real, 'r-', c= aotomat_green,
         label=r'Re($\varepsilon$)')
plt.plot(lambdas , epsilon.imag, 'r--', c= aotomat_purple,
         label=r'Im($\varepsilon$)')
plt.xlabel(r'$\lambda$ ($\mu m$)')
name = yamlFile[5:][:-4]
plt.title("complex permittivity from " + name)
plt.legend(loc=0)
plt.show()


############################################################################
# .. [JC1972]
#          (P. B. Johnson and R. W. Christy. Optical constants of the noble metals, Phys. Rev. B 6, 4370-4379 (1972)).
# .. _`Johnson and Christy`:
#     https://doi.org/10.1103/PhysRevB.6.4370
