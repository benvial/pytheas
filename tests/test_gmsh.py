
import os


def test_gmsh():
    rc = os.system("gmsh --info")
    assert rc == 0
