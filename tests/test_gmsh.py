import os


def test_gmsh():
    print("\n")
    rc = os.system("gmsh --info")
    proc = subprocess.Popen(["getdp --info"], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    assert err is None


def test_t1():
    rc = os.system("gmsh t1.geo -2 -v 8")
    assert rc == 0
    os.system("rm t1.msh")
