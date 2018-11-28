#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


import pygmsh
import meshio


geom = pygmsh.built_in.Geometry()

lcar = 0.1
h = 25
w = 10
length = 100
# x_fin = -0.5 * length
cr = 1

f = 0.5 * w
y = [-f, -f + cr, +f - cr, +f]
z = [0.0, 0, 0]
f = 0.5 * cr
x = [-f, f]
points = []
points.append(geom.add_point((x[0], y[0], z[0]), lcar=lcar))
points.append(geom.add_point((x[1], y[0], z[0]), lcar=lcar))

lines = []
lines.append(geom.add_line(points[0], points[1]))

line, surf, _ = geom.extrude(lines[0], translation_axis=[0,1,0], num_layers=7,recombine=True)
# line_loop = geom.add_line_loop(l1[1])

geom.extrude(surf, translation_axis=[0, 0, 1], num_layers=7,recombine=True)

# geom = pygmsh.opencascade.Geometry()

# 
# line_loop = geom.add_line_loop(lines)
# surface = geom.add_plane_surface(line_loop)
# geom.extrude(surface, translation_axis=[length, 0, 0])

# geom.add_physical_volume(bkg)
geom.add_raw_code('Coherence;')

code_built_in = geom.get_code().replace("\'", "")

geom = pygmsh.opencascade.Geometry(    )

des = geom.add_rectangle([1,1,1], 1,2)
_,des = geom.extrude(des, [0,0,1])

code_opencascade = geom.get_code().replace("\'", "")


code = code_built_in + "\n" + code_opencascade

geo_filename="ext.geo"

with open(geo_filename, "w") as f:
    f.write(code)
