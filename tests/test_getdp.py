# test for getdp
import subprocess
from pytheas.tools.femio import getdp


def test_getdp():
    out = subprocess.call([getdp, "--info"])
    assert out == 0
