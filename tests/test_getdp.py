# test for getdp
import subprocess


def test_getdp():
    print("\n")
    proc = subprocess.Popen(["getdp --info"], stdout=subprocess.PIPE, shell=True)
    _, err = proc.communicate()
    assert err is None
