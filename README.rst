
.. image:: https://travis-ci.org/benvial/pytheas.svg?branch=master
   :target: https://travis-ci.org/benvial/pytheas
   :alt: Travis CI build status (Linux)

.. image:: https://codecov.io/github/benvial/pytheas/coverage.svg?branch=master
   :target: https://codecov.io/github/benvial/pytheas?branch=master
   :alt: Code Coverage

Python Electromagnetic Analysis and Simulation with the Finite Element Method
-----------------------------------------------------------------------------

.. inclusion-marker-do-not-remove

Pytheas is a `Python <http://www.python.org/>`_ package for creating,
running and postprocessing electrodynamic simulations. It is based on open
source software `Gmsh <http://www.gmsh.info/>`_. for creating
geometries and mesh generation, and `GetDP <http://www.getdp.info/>`_. for solving
the underlying partial differential equations with the finite
element method.

It features built in models of:

- periodic media in 2D and 3D with computation of diffraction efficiencies
- scattering analysis in 2D and 3D
- bloch mode analysis of metamaterials
- treatment of open geometries with perfectly matched layers
- tools to define arbitrary permittivity distributions
- quasi-normal mode analysis
- two scale convergence homogenization
- tools for topology optimization in 2D
- built-in refractive index database
