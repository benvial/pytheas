"""
Tools for gmsh/getdp control and input/output.

"""


import meshio
import os
import subprocess
import numpy as np
import shutil
from pyonelab import gmsh, getdp


def make_var_str(varname, varvalue):
    s = varname + " = " + str(varvalue) + ";"
    return s


def append_inputs(str_to_append, data_file):
    f = open(data_file, "a")
    f.write("\n")
    f.write(str_to_append)
    f.close()


def make_inputs(data_dict):
    S = ""
    for key in data_dict:
        s = "\n" + make_var_str(key, data_dict[key])
        S += s
    return S


def write_inputs(string, data_file, close=True):
    f = open(data_file, "w")
    f.write(string)
    if close:
        f.close()
    return f


def get_content(filename):
    with open(filename) as f:
        conts = f.read()
        f.close()
    return conts


def maketmp(content, filename, dirname="", mode="w"):
    path = os.path.join(dirname, filename)
    with open(path, mode) as f:
        f.write(content)
        f.close()
    return path


def mesh_model(
    path_mesh, path_geo, mesh_format="msh2", dim=None, verbose=0, other_option=""
):
    """Mesh the model using Gmsh_

    .. _Gmsh:
        http://gmsh.info/
    """
    dim = dim or [1, 2]
    list_dim = ["-{}".format(d) for d in dim]
    command = [gmsh, "-v", str(verbose), "-format", mesh_format, path_geo]
    command += list_dim + other_option.split() + ["-o", path_mesh]
    subprocess.call(command)


def solve_problem(resolution, path_pro, path_mesh, path_pos=None, verbose=0, argstr=""):
    command = [
        getdp,
        "-v",
        str(verbose),
        path_pro,
        "-pre",
        resolution,
        "-msh",
        path_mesh,
    ]
    if path_pos:
        command += ["-gmshread"] + path_pos.split()
    command += ["-cal", "-v2"] + argstr.split()
    subprocess.call(command)


def make_content_mesh_pos(nodes, els, dom, celltype):
    nodes_ID, nodes_coords = nodes
    els_ID, _, els_nodes_ID, geom_ID_dom = els
    s = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
    nnodes = len(nodes_ID)
    s += "$Nodes\n"
    s += str(nnodes) + "\n"
    for i in range(nnodes):
        x, y, z = nodes_coords[i]
        s += str(nodes_ID[i]) + " " + str(x) + " " + str(y) + " " + str(z) + "\n"
    s += "$EndNodes\n"
    nels = len(els_ID)
    s += "$Elements\n"
    s += str(nels) + "\n"
    for i in range(nels):
        if celltype is "triangle":
            s1 = str(2)
        elif celltype is "quad":
            s1 = str(3)
        elif celltype is "tetra":
            s1 = str(4)
        s += (
            str(els_ID[i])
            + " "
            + s1
            + " 2 "
            + str(dom)
            + " "
            + str(geom_ID_dom[i])
            + " "
        )
        for j in els_nodes_ID[i]:
            s += str(j) + " "
        s += "\n"
    s += "$EndElements\n"
    return s


def get_nodes(path_mesh, physical_ID, celltype):
    mesh = meshio.read(path_mesh)
    # phys_ID = mesh.cell_data[celltype]["gmsh:physical"]
    types_cells = [_.type for _ in mesh.cells]
    celltype_index = [i for i, _ in enumerate(types_cells) if _ == celltype][0]
    phys_ID = mesh.cell_data["gmsh:physical"][celltype_index]
    domain = phys_ID == physical_ID
    cell = mesh.cells[celltype_index]
    els_nodes_ID = cell.data[domain]
    nodes_ID_domain = np.unique(els_nodes_ID.flatten())
    nodes_coords_domain = mesh.points[nodes_ID_domain]
    return nodes_ID_domain + 1, nodes_coords_domain


def get_elements(path_mesh, physical_ID, celltype):
    mesh = meshio.read(path_mesh)
    types_cells = [_.type for _ in mesh.cells]
    celltype_index = [i for i, _ in enumerate(types_cells) if _ == celltype][0]
    phys_ID = mesh.cell_data["gmsh:physical"][celltype_index]
    domain = phys_ID == physical_ID
    n = 1
    for ik, cell in enumerate(mesh.cells):
        if cell.type is celltype:
            n += 0
        else:
            n += len(mesh.cell_data["gmsh:physical"][ik])
    el_ID = np.arange(0, len(domain))[domain] + n
    cell = mesh.cells[celltype_index].data
    els_nodes_ID = cell[domain]
    el_center = np.mean(mesh.points[els_nodes_ID], axis=1)
    return el_ID, el_center, els_nodes_ID + 1, None


def make_pos(ID, data, content_mesh, viewname, celltype="nodes", mesh_format=2):
    s = content_mesh
    if not s:
        if mesh_format == 2:
            meshversion = "2.2 0 8"
        else:
            meshversion = "4.1 0 8"
        s = "$MeshFormat\n{}\n$EndMeshFormat\n".format(meshversion)
    N = len(ID)
    if celltype is "nodes":
        str_start, str_end = "$NodeData\n", "$EndNodeData\n"
    elif celltype is "elements":
        str_start, str_end = "$ElementData\n", "$EndElementData\n"
    elif celltype is "elements_nodes":
        str_start, str_end = "$ElementNodeData\n", "$EndElementNodeData\n"
    else:
        raise TypeError("Wrong celltype specified: choose between nodes and elements")
    for dat, name in zip([data.real, data.imag], ["_real", "_imag"]):
        s += str_start
        s += '1\n"' + viewname + name + '"\n1\n0\n3\n0\n1\n' + str(N) + "\n"
        for idf, value in zip(ID, dat):
            if celltype is "elements_nodes":
                s += str(int(idf)) + " 3 " + (str(value.real) + " ") * 3 + "\n"
            else:
                s += str(int(idf)) + " " + str(value.real) + "\n"
        s += str_end
    return s


def open_gmsh(path_mesh, path_geo, pos_list=None, verbose=2):
    pos_list = pos_list or []
    command = [gmsh, path_geo, path_mesh] + pos_list + ["-n", "-v", str(verbose), "&"]

    # subprocess.call(command, shell=True)
    subprocess.call(" ".join(command), shell=True)
    # os.system(" ".join(command))


def postpro_commands(postop, path_pro, path_mesh, path_pos=None, verbose=0):
    """Generate a command list for postprocessing by GetDP (see main.pro
    file in ./base folder for default available postprocessings, or to add
    your own)

    Args:
        postop (str): The name of the postoperation to perform.
        path_pro (str): Path to the .pro file
        path_mesh (str): Path to the .msh file
        path_pos (str , optional): Path to a file to be read by ``gmshread``.
        verbose (int): verbosity level
        Defaults to None.

    Returns:
        list : The list of strings to be oscommanded.


    """
    path_res = path_pro[:-4] + ".res"
    s = (
        getdp
        + " -v "
        + str(verbose)
        + " "
        + path_pro
        + " -res "
        + path_res
        + " -msh "
        + path_mesh
    )
    if path_pos:
        s += " -gmshread " + path_pos
    s += " -pos " + postop + " -v2"
    return s.split()


def load_node_table(filename):
    nodenumber = np.loadtxt(filename, usecols=[0], skiprows=1)
    values = np.loadtxt(filename, usecols=[1], skiprows=1) + 1j * np.loadtxt(
        filename, usecols=[2], skiprows=1
    )
    return nodenumber, values


def load_table(filename):
    return np.loadtxt(filename, usecols=[3]) + 1j * np.loadtxt(filename, usecols=[4])


# def load_node_table_vect(filename):
#     nodenumber = np.loadtxt(filename, usecols=[0], skiprows=1)
#     values = np.loadtxt(filename, usecols=[1], skiprows=1) + 1j * np.loadtxt(
#         filename, usecols=[2], skiprows=1
#     )
#     return nodenumber, values


def load_table_vect(filename):
    vect = []
    for i in range(3):
        comp = 0
        for j in range(2):
            comp += np.loadtxt(filename, usecols=[3 + i + j * 3]) * np.exp(
                j * 1j * np.pi / 2
            )
        vect.append(comp)
    return vect


def load_table_tens(filename):
    a = np.loadtxt(filename)[3:]
    return (a[:9] + 1j * a[9:]).reshape((3, 3))


def load_node_table_vect(filename):
    nodenumber = np.loadtxt(filename, usecols=[0], skiprows=1)
    vect = []
    for i in range(3):
        comp = 0
        for j in range(2):
            comp += np.loadtxt(filename, usecols=[1 + i + j * 3], skiprows=1) * np.exp(
                j * 1j * np.pi / 2
            )
        vect.append(comp)
    return nodenumber, vect


def load_element_table_vect(filename):
    # nodenumber = np.loadtxt(filename, usecols=[0], skiprows=1)
    vect = []
    for i in range(3):
        comp = 0
        for j in range(2):
            comp += np.loadtxt(filename, usecols=[1 + i + j * 3], skiprows=1) * np.exp(
                j * 1j * np.pi / 2
            )
        vect.append(comp)
    # return nodenumber, vect
    return vect


def load_timetable(filename):
    return np.loadtxt(filename, usecols=[5]) + 1j * np.loadtxt(filename, usecols=[6])


def load_ev_timetable(filename):
    return np.loadtxt(filename, usecols=[1]) + 1j * np.loadtxt(filename, usecols=[5])


def points2geo(points, lc_incl, output_path="./tmp.geo", startpoint=1000):
    x, y = points
    fout = open(output_path, "w")
    # lc_name = "%s_lc" % filename[0:3]
    # Format
    # Point(1) = {0, 0, 0, lc};
    # fout.write("%s = 0.005;\n" % lc_name)
    j = startpoint
    n_lines = len(x)
    for i in range(n_lines):
        outputline = "Point(%i) = { %8.8f, %8.8f, %8.8f, %s};\n " % (
            j,
            x[i],
            y[i],
            0,
            lc_incl,
        )
        j = j + 1
        fout.write(outputline)
    # gmsh bspline format
    # Write out splinefit line
    fout.write(
        "BSpline(%i) = {%i:%i};\n" % (startpoint, startpoint, startpoint + n_lines - 1)
    )
    fout.close()
