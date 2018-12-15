#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT

"""
The :mod:`pytheas.periodic2D` module implements the resolution of the
scalar wave equation for TE and TM polarization for
mono-periodic stuctures in 2D:

- subject to an incident plane wave (diffraction problem) with calculation
  of the diffraction efficiencies, absorption and energy balance.
- eigenvalues and eigenmodes (modal analysis)

"""


from .femmodel import Periodic2D
