import os
import subprocess
from pytheas.tools.femio import gmsh


def test_gmsh():
    out = subprocess.call([gmsh, "--info"])
    assert out == 0


def test_t1():
    out = subprocess.call([gmsh, "t1.geo", "-2"])
    print(out)
    assert out == 0
    assert os.path.isfile("t1.msh")
    os.remove("t1.msh")
