#!/usr/bin/env python

import numpy as np

pi = np.pi


def ellipse(Rinclx, Rincly, x0, y0, rot_incl, nt=360):
    c, s = np.cos(rot_incl), np.sin(rot_incl)
    Rot = np.array([[c, -s], [s, c]])
    theta = np.linspace(-pi, pi, nt)
    x = Rinclx * np.sin(theta)
    y = Rincly * np.cos(theta)
    x, y = np.linalg.linalg.dot(Rot, np.array([x, y]))
    points = x + x0, y + y0
    return points


def circle(R, x0, y0, nt=360):
    return ellipse(R, R, x0, y0, 0, nt=nt)
