import nlopt
import sys
import numpy as np
import numpy.testing as npt


def test_nlopt_import():
    assert "nlopt" in sys.modules


def myfunc(x, grad):
    if grad.size > 0:
        grad[0] = 0.0
        grad[1] = 0.5 / np.sqrt(x[1])
    return np.sqrt(x[1])


def myconstraint(x, grad, a, b):
    if grad.size > 0:
        grad[0] = 3 * a * (a * x[0] + b) ** 2
        grad[1] = -1.0
    return (a * x[0] + b) ** 3 - x[1]


def test_nlopt():
    opt = nlopt.opt(nlopt.LD_MMA, 2)
    opt.set_lower_bounds([-float("inf"), 0])
    opt.set_min_objective(myfunc)
    opt.add_inequality_constraint(lambda x, grad: myconstraint(x, grad, 2, 0), 1e-8)
    opt.add_inequality_constraint(lambda x, grad: myconstraint(x, grad, -1, 1), 1e-8)
    opt.set_xtol_rel(1e-4)
    x = opt.optimize([1.234, 5.678])
    minf = opt.last_optimum_value()
    # numevals = opt.get_numevals()
    res = opt.last_optimize_result()
    print("optimum at ", x[0], x[1])
    print("minimum value = ", minf)
    print("result code = ", res)
    # print("nevals = ", numevals)
    min_fref = 0.5443310476200902
    xref = np.array([0.3333333346933468, 0.29629628940318486])
    assert res == 4
    # assert numevals == 11
    assert minf == min_fref
    npt.assert_almost_equal(xref, x, decimal=3)
