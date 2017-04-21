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


def assemble_Hk(bfgs_in):
    """
    Assemble Hk
    """
    Hy, e = dl.Vector(), dl.Vector()
    bfgs_in.model.Prior.init_vector(e, 0)
    bfgs_in.model.Prior.init_vector(Hy, 0)
    dim = len(e.array())
    Hk = []
    for ii in range(dim):
        e.zero()
        e[ii] = 1.0
        bfgs_in.apply_Hk(e, Hy)
        Hk.append(Hy.array())
    return np.array(Hk)

def test2():
    """
    Test that assembled matrix Hk is SPD
    """
    mesh = dl.UnitSquareMesh(10, 10)
    V = dl.FunctionSpace(mesh, 'CG', 1)
    zp = ZeroPrior(V)
    model = Model([], zp, [])

    s, y, = dl.Vector(), dl.Vector()
    zp.init_vector(s, 0)
    zp.init_vector(y, 0)
    dim = len(s.array())
    print 'dim={}'.format(dim)

    for ii in range(10):
        print '\nTest {}'.format(ii)
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

        Hk = assemble_Hk(bfgs)

        diffHkHkT = np.linalg.norm(Hk - Hk.T)/np.linalg.norm(Hk)
        print 'Error due to asymetry of Hk = {}'.format(diffHkHkT),
        if diffHkHkT < 1e-12:   print '\t==>> OK!!'
        else:   print '\t\t\t<<<<==== PROBLEM!?'

        eig = np.linalg.eigvalsh(Hk)
        print 'smallest and largest eigenvalues = {}, {}'.format(np.min(eig), np.max(eig)),
        if np.min(eig) > 0.0:   print '\t==>> OK!!'
        else:   print '\t\t\t<<<<==== PROBLEM!?'



def test3():
    """
    Compare with assembling matrix Hk directly
    """
    mesh = dl.UnitSquareMesh(10, 10)
    V = dl.FunctionSpace(mesh, 'CG', 1)
    zp = ZeroPrior(V)
    model = Model([], zp, [])

    s, y, = dl.Vector(), dl.Vector()
    zp.init_vector(s, 0)
    zp.init_vector(y, 0)
    dim = len(s.array())
    print 'dim={}'.format(dim)

    for ii in range(10):
        print '\nTest {}'.format(ii)
        bfgs = BFGS(model)
        Hm = np.eye(dim)

        for ii in range(30):
            s.set_local((ii+1.0)*np.random.randn(dim))
            y.set_local(np.random.randn(dim))
            # for BFGS, we need to have s^T.y > 0
            if s.inner(y) < 0.0:
                y *= -1.0
            bfgs.S.append(s.copy())
            bfgs.Y.append(y.copy())
            rho = 1./s.inner(y)
            bfgs.R.append(rho)

            # Update matrix form
            sk = s.array().reshape((dim,-1))
            yk = y.array().reshape((dim,-1))
            Ir = np.eye(dim) - rho*sk.dot(yk.T)
            Hm = Ir.dot(Hm.dot(Ir.T)) + rho*sk.dot(sk.T)

        Hk = assemble_Hk(bfgs)

        diffHkHm = np.linalg.norm(Hk-Hm)/np.linalg.norm(Hk)
        print 'Diff between Hk and Hm = {}'.format(diffHkHm),
        if diffHkHm < 1e-12:   print '\t==>> OK!!'
        else:   print '\t\t\t<<<<==== PROBLEM!?'



#TODO:def test4():
    """
    Check how quickly BFGS Hessian converges toward Hessian of quadratic form
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
    elif case == 3:
        test3()
