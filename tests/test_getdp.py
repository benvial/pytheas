# test for getdp
import subprocess


def test_getdp():
    proc = subprocess.Popen(["getdp", "--info"], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    print("output: \n", out)
    print("error: \n", err)
    assert err is None
