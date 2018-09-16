
.. image:: https://img.shields.io/travis/benvial/pytheas/master.svg?style=for-the-badge
   :target: https://travis-ci.org/benvial/pytheas
   :alt: Travis CI build status (Linux)

.. image:: https://img.shields.io/codecov/c/github/benvial/pytheas.svg?style=for-the-badge
   :target: https://codecov.io/github/benvial/pytheas?branch=master
   :alt: Code Coverage

.. image:: https://img.shields.io/codacy/grade/e27821fb6289410b8f58338c7e0bc686.svg?style=for-the-badge
   :target: https://app.codacy.com/app/benvial/pytheas/dashboard
   :alt: Codacy grade

.. image:: https://img.shields.io/github/license/mashape/apistatus.svg?style=for-the-badge
   :alt: Licence: MIT

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg?style=for-the-badge
   :alt: Code style: black
   

Python Electromagnetic Analysis and Simulation with the Finite Element Method
-----------------------------------------------------------------------------

.. inclusion-marker-do-not-remove

Pytheas is a `Python <http://www.python.org/>`_ package for creating,
running and postprocessing electrodynamic simulations. It is based on open
source software `Gmsh <http://www.gmsh.info/>`_ for creating
geometries and mesh generation, and `GetDP <http://www.getdp.info/>`_ for solving
the underlying partial differential equations with the finite
element method.

It features built in models of:

- periodic media in 2D and 3D with computation of diffraction efficiencies
- scattering analysis in 2D and 3D
- Bloch mode analysis of metamaterials
- treatment of open geometries with perfectly matched layers
- tools to define arbitrary permittivity distributions
- quasi-normal mode analysis
- two scale convergence homogenization
- tools for topology optimization in 2D
- built-in refractive index database
