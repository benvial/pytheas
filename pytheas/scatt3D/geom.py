#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


import pygmsh
import os


def make_geom(ext=False):
    pml_c = [None for _ in range(8)]
    pml_x = [None for _ in range(2)]
    pml_y = [None for _ in range(2)]
    pml_z = [None for _ in range(2)]
    pml_xy = [None for _ in range(4)]
    pml_xz = [None for _ in range(4)]
    pml_yz = [None for _ in range(4)]

    geom_oc = pygmsh.opencascade.Geometry()

    geom_oc.add_raw_code('Include "parameters.dat";')
    geom_oc.add_raw_code("lc_des = lambda0/(Sqrt[eps_des_re]*parmesh_des);")
    geom_oc.add_raw_code("lc_host = lambda0/(Sqrt[eps_host_re]*parmesh);")
    geom_oc.add_raw_code("lc_pml = lambda0/(Sqrt[eps_host_re]*parmesh_pml);")
    geom_oc.add_raw_code("lc_sph = lambda0/(Sqrt[eps_sph_re]*parmesh_incl);")

    des = geom_oc.add_rectangle(
        ["-hx_des/2", "-hy_des/2", "-hz_des/2"],
        "hx_des",
        "hy_des",
        char_length="lc_des",
    )
    _, des, _ = geom_oc.extrude(des, [0, 0, "hz_des"])

    box = geom_oc.add_rectangle(
        ["-hx_box/2", "-hy_box/2", "-hz_box/2"],
        "hx_box",
        "hy_box",
        char_length="lc_host",
    )
    _, box, _ = geom_oc.extrude(box, [0, 0, "hz_box"])
    # union = geom_oc.boolean_difference([box_], [des])

    pml1 = geom_oc.add_rectangle(
        ["-h_pml-hx_box/2", "-h_pml-hy_box/2", "-h_pml-hz_box/2"],
        "h_pml",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_c[0], _ = geom_oc.extrude(pml1, [0, 0, "h_pml"])
    pml2 = geom_oc.add_rectangle(
        ["-h_pml-hx_box/2", "-h_pml-hy_box/2", "-hz_box/2"],
        "h_pml",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_xy[0], _ = geom_oc.extrude(pml2, [0, 0, "hz_box"])
    pml3 = geom_oc.add_rectangle(
        ["-h_pml-hx_box/2", "-h_pml-hy_box/2", "hz_box/2"],
        "h_pml",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_c[1], _ = geom_oc.extrude(pml3, [0, 0, "h_pml"])

    pml1 = geom_oc.add_rectangle(
        ["hx_box/2", "-h_pml-hy_box/2", "-h_pml-hz_box/2"],
        "h_pml",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_c[2], _ = geom_oc.extrude(pml1, [0, 0, "h_pml"])
    pml2 = geom_oc.add_rectangle(
        ["hx_box/2", "-h_pml-hy_box/2", "-hz_box/2"],
        "h_pml",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_xy[1], _ = geom_oc.extrude(pml2, [0, 0, "hz_box"])
    pml3 = geom_oc.add_rectangle(
        ["hx_box/2", "-h_pml-hy_box/2", "hz_box/2"],
        "h_pml",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_c[3], _ = geom_oc.extrude(pml3, [0, 0, "h_pml"])

    pml1 = geom_oc.add_rectangle(
        ["-h_pml-hx_box/2", "hy_box/2", "-h_pml-hz_box/2"],
        "h_pml",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_c[4], _ = geom_oc.extrude(pml1, [0, 0, "h_pml"])
    pml2 = geom_oc.add_rectangle(
        ["-h_pml-hx_box/2", "hy_box/2", "-hz_box/2"],
        "h_pml",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_xy[2], _ = geom_oc.extrude(pml2, [0, 0, "hz_box"])
    pml3 = geom_oc.add_rectangle(
        ["-h_pml-hx_box/2", "hy_box/2", "hz_box/2"],
        "h_pml",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_c[5], _ = geom_oc.extrude(pml3, [0, 0, "h_pml"])

    pml1 = geom_oc.add_rectangle(
        ["hx_box/2", "hy_box/2", "-h_pml-hz_box/2"],
        "h_pml",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_c[6], _ = geom_oc.extrude(pml1, [0, 0, "h_pml"])
    pml2 = geom_oc.add_rectangle(
        ["hx_box/2", "hy_box/2", "-hz_box/2"], "h_pml", "h_pml", char_length="lc_pml"
    )
    _, pml_xy[3], _ = geom_oc.extrude(pml2, [0, 0, "hz_box"])
    pml3 = geom_oc.add_rectangle(
        ["hx_box/2", "hy_box/2", "hz_box/2"], "h_pml", "h_pml", char_length="lc_pml"
    )
    _, pml_c[7], _ = geom_oc.extrude(pml3, [0, 0, "h_pml"])

    pml1 = geom_oc.add_rectangle(
        ["-hx_box/2", "-h_pml-hy_box/2", "-h_pml-hz_box/2"],
        "hx_box",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_yz[0], _ = geom_oc.extrude(pml1, [0, 0, "h_pml"])
    pml2 = geom_oc.add_rectangle(
        ["-hx_box/2", "-h_pml-hy_box/2", "-hz_box/2"],
        "hx_box",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_y[0], _ = geom_oc.extrude(pml2, [0, 0, "hz_box"])
    pml3 = geom_oc.add_rectangle(
        ["-hx_box/2", "-h_pml-hy_box/2", "hz_box/2"],
        "hx_box",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_yz[1], _ = geom_oc.extrude(pml3, [0, 0, "h_pml"])

    pml1 = geom_oc.add_rectangle(
        ["-hx_box/2", "hy_box/2", "-h_pml-hz_box/2"],
        "hx_box",
        "h_pml",
        char_length="lc_pml",
    )
    _, pml_yz[2], _ = geom_oc.extrude(pml1, [0, 0, "h_pml"])
    pml2 = geom_oc.add_rectangle(
        ["-hx_box/2", "hy_box/2", "-hz_box/2"], "hx_box", "h_pml", char_length="lc_pml"
    )
    _, pml_y[1], _ = geom_oc.extrude(pml2, [0, 0, "hz_box"])
    pml3 = geom_oc.add_rectangle(
        ["-hx_box/2", "hy_box/2", "hz_box/2"], "hx_box", "h_pml", char_length="lc_pml"
    )
    _, pml_yz[3], _ = geom_oc.extrude(pml3, [0, 0, "h_pml"])

    pml1 = geom_oc.add_rectangle(
        ["-h_pml-hx_box/2", "-hy_box/2", "-h_pml-hz_box/2"],
        "h_pml",
        "hy_box",
        char_length="lc_pml",
    )
    _, pml_xz[0], _ = geom_oc.extrude(pml1, [0, 0, "h_pml"])
    pml2 = geom_oc.add_rectangle(
        ["-h_pml-hx_box/2", "-hy_box/2", "-hz_box/2"],
        "h_pml",
        "hy_box",
        char_length="lc_pml",
    )
    _, pml_x[0], _ = geom_oc.extrude(pml2, [0, 0, "hz_box"])
    pml3 = geom_oc.add_rectangle(
        ["-h_pml-hx_box/2", "-hy_box/2", "hz_box/2"],
        "h_pml",
        "hy_box",
        char_length="lc_pml",
    )
    _, pml_xz[1], _ = geom_oc.extrude(pml3, [0, 0, "h_pml"])

    pml1 = geom_oc.add_rectangle(
        ["hx_box/2", "-hy_box/2", "-h_pml-hz_box/2"],
        "h_pml",
        "hy_box",
        char_length="lc_pml",
    )
    _, pml_xz[2], _ = geom_oc.extrude(pml1, [0, 0, "h_pml"])
    pml2 = geom_oc.add_rectangle(
        ["hx_box/2", "-hy_box/2", "-hz_box/2"], "h_pml", "hy_box", char_length="lc_pml"
    )
    _, pml_x[1], _ = geom_oc.extrude(pml2, [0, 0, "hz_box"])
    pml3 = geom_oc.add_rectangle(
        ["hx_box/2", "-hy_box/2", "hz_box/2"], "h_pml", "hy_box", char_length="lc_pml"
    )
    _, pml_xz[3], _ = geom_oc.extrude(pml3, [0, 0, "h_pml"])
    #

    pml1 = geom_oc.add_rectangle(
        ["-hx_box/2", "-hy_box/2", "-h_pml-hz_box/2"],
        "hx_box",
        "hy_box",
        char_length="lc_pml",
    )
    _, pml_z[0], _ = geom_oc.extrude(pml1, [0, 0, "h_pml"])
    pml2 = geom_oc.add_rectangle(
        ["-hx_box/2", "-hy_box/2", "hz_box/2"], "hx_box", "hy_box", char_length="lc_pml"
    )
    _, pml_z[1], _ = geom_oc.extrude(pml2, [0, 0, "h_pml"])

    sph = geom_oc.add_ball(["x_sph", "y_sph", "z_sph"], "R_sph", char_length="lc_sph")
    bkg = geom_oc.boolean_difference([box], [sph, des])
    sph = geom_oc.add_ball(["x_sph", "y_sph", "z_sph"], "R_sph", char_length="lc_sph")
    #

    if ext:
        #
        # geom = pygmsh.built_in.Geometry()
        # points = []
        # points.append(geom.add_point(("-hx_des/2", "-hy_des/2","-hz_des/2"), lcar="lc_des"))
        # points.append(geom.add_point(("hx_des/2", "-hy_des/2","-hz_des/2"), lcar="lc_des"))
        # line = geom.add_line(points[0], points[1])
        # lines, surface, _ = geom.extrude(line, translation_axis=[0,"hy_des",0], num_layers=7,recombine=True)
        #
        # # line_loop = geom.add_line_loop(lines)
        # # surface = geom.add_plane_surface(line_loop)
        #
        # _,des,_ = geom.extrude(surface, translation_axis=[0, 0, "hz_des"], num_layers=7,recombine=True)
        #
        # code= geom.get_code().replace("\'", "")
        # print(code)
        #
        geom_oc.add_raw_code('SetFactory("Built-In");')
        geom_oc.add_raw_code(
            "p135 = newp; Point(p135) = {-hx_des/2, -hy_des/2, -hz_des/2, lc_des};"
        )
        geom_oc.add_raw_code(
            "p136 = newp; Point(p136) = {hx_des/2, -hy_des/2, -hz_des/2, lc_des};"
        )
        geom_oc.add_raw_code("l46 = newl; Line(l46) = {p135, p136};")

        geom_oc.add_raw_code(
            "ex1[] = Extrude {0,hy_des,0} {Line{l46}; Layers{Ceil(hy_des/lc_des)}; Recombine;};"
        )
        geom_oc.add_raw_code(
            "des_dom[] = Extrude {0,0,hz_des} {Surface{ex1[1]}; Layers{Ceil(hz_des/lc_des)}; Recombine;};"
        )
        des.id = "des_dom[1]"
    else:
        des = geom_oc.add_rectangle(
            ["-hx_des/2", "-hy_des/2", "-hz_des/2"],
            "hx_des",
            "hy_des",
            char_length="lc_des",
        )
        _, des, _ = geom_oc.extrude(des, [0, 0, "hz_des"])

    printpoint = geom_oc.add_point([0, 0, 0])

    geom_oc.add_physical(bkg)
    geom_oc.add_physical(des)
    geom_oc.add_physical(sph)
    geom_oc.add_physical(pml_c)
    geom_oc.add_physical(pml_xy)
    geom_oc.add_physical(pml_yz)
    geom_oc.add_physical(pml_xz)
    geom_oc.add_physical(pml_x)
    geom_oc.add_physical(pml_y)
    geom_oc.add_physical(pml_z)
    geom_oc.add_physical(printpoint)

    geom_oc.add_raw_code("Coherence;")
    geom_oc.add_raw_code("Coherence;")
    geom_oc.add_raw_code("Coherence;")

    geom_oc.add_raw_code("Mesh.Algorithm=6;")
    geom_oc.add_raw_code("Mesh.Algorithm3D=4;")

    code = geom_oc.get_code().replace("'", "")

    dir_path = os.path.dirname(os.path.abspath(__file__))

    geo_filename = dir_path + "/base/geometry.geo"
    with open(geo_filename, "w") as f:
        f.write(code)


# import os
#
# os.system("gmsh " + geo_filename)
