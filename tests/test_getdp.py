
import os


def test_getdp():
    rc = os.system("getdp --info")
    assert rc == 0
