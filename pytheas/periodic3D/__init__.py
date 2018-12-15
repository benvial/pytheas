#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT

"""
The :mod:`pytheas.periodic3D` module implements the resolution of the
vectorial wave equation in 3D for bi-periodic structure:

- subject to an incident plane wave (scattering problem), with calculation
  of the diffraction efficiencies, absorption and energy balance.
- eigenvalues and eigenmodes (modal analysis)

"""


from .femmodel import Periodic3D
