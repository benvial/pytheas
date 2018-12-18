#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT

"""
The :mod:`pytheas.scatt3D` module implements the resolution of the
vectorial wave equation in 3D:

- subject to an incident plane wave (scattering problem)
- eigenvalues and eigenmodes (modal analysis)

"""


from .femmodel import Scatt3D
