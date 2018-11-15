import os


def test_gmsh():
    rc = os.system("gmsh --info")
    assert rc == 0


def test_t1():
    rc = os.system("gmsh t1.geo -2 -v 8")
    assert rc == 0
    os.system("rm t1.msh")
