# test for getdp
import subprocess


def test_getdp():
    out = subprocess.call(["getdp", "--info"])
    assert out == 0
