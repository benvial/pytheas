import os
import subprocess


def test_gmsh():
    out = subprocess.call(["gmsh", "--info"])
    assert out == 0


def test_t1():
    out = subprocess.call(["gmsh", "t1.geo", "-2"])
    print(out)
    subprocess.call(["rm", "t1.msh"])
    assert out == 0
