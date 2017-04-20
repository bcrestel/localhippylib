"""
Test application of the matrix H in BFGS
"""

import sys
import numpy as np
import dolfin as dl
from hippylib import *

def test1():
    """
    Test that assembled matrix Hk verify: Hk+1 * yk = sk
    """
    mesh = dl.UnitSquareMesh(10, 10)
    V = dl.FunctionSpace(mesh, 'CG', 1)
    zp = ZeroPrior(V)
    model = Model([], zp, [])

    s, y, Hy = dl.Vector(), dl.Vector(), dl.Vector()
    zp.init_vector(s, 0)
    zp.init_vector(y, 0)
    zp.init_vector(Hy, 0)
    dim = len(s.array())

    bfgs = BFGS(model)

    for ii in range(10):
        s.set_local((ii+1.0)*np.random.randn(dim))
        y.set_local(np.random.randn(dim))
        # for BFGS, we need to have s^T.y > 0
        if s.inner(y) < 0.0:
            y *= -1.0
        bfgs.S.append(s.copy())
        bfgs.Y.append(y.copy())
        bfgs.R.append(1./s.inner(y))

        bfgs.apply_Hk(y, Hy)
        ns = dl.norm(s)
        ny = dl.norm(y)
        nHy = dl.norm(Hy)
        err = dl.norm(Hy-s)/ns
        print '|s|={}, |Hy|={}, |y|={}, |Hy-s|/|s|={:.2e}'.format(\
        ns, nHy, ny, err),
        if err < 1e-12:
            print '\t==>> OK!!'



def test2():
    """
    Test that assembled matrix Hk is SPD
    """
    mesh = dl.UnitSquareMesh(10, 10)
    V = dl.FunctionSpace(mesh, 'CG', 1)
    zp = ZeroPrior(V)
    model = Model([], zp, [])

    s, y, Hy, e = dl.Vector(), dl.Vector(), dl.Vector(), dl.Vector()
    zp.init_vector(s, 0)
    zp.init_vector(y, 0)
    zp.init_vector(Hy, 0)
    zp.init_vector(e, 0)
    dim = len(s.array())

    bfgs = BFGS(model)

    for ii in range(10):
        s.set_local((ii+1.0)*np.random.randn(dim))
        y.set_local(np.random.randn(dim))
        # for BFGS, we need to have s^T.y > 0
        if s.inner(y) < 0.0:
            y *= -1.0
        bfgs.S.append(s.copy())
        bfgs.Y.append(y.copy())
        bfgs.R.append(1./s.inner(y))

    # Assemble Hk
    Hk = []
    for ii in range(dim):
        e.zero()
        e[ii] = 1.0
        bfgs.apply_Hk(e, Hy)
        Hk.append(Hy.array())
    Hk = np.array(Hk)

    diffHkHkT = np.linalg.norm(Hk - Hk.T)/np.linalg.norm(Hk)
    print 'Error due to asymetry of Hk = {}'.format(diffHkHkT)

    eig = np.linalg.eigvalsh(Hk)
    print 'smallest and largest eigenvalues = {}, {}'.format(np.min(eig), np.max(eig))



#TODO
#def test3():
    """
    Test BFGS on simple convex optimization problem
    """



if __name__ == "__main__":
    try:
        case = int(sys.argv[1])
    except:
        print 'Usage:\n\tpython {} <test_case_nb>'.format(sys.argv[0])
        sys.exit(1)

    if case == 1:
        test1()
    elif case == 2:
        test2()
