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

    bfgs = BFGS_operator()

    for ii in range(10):
        s.set_local((ii+1.0)*np.random.randn(dim))
        y.set_local(np.random.randn(dim))
        bfgs.update(s, y)

        bfgs.solve(Hy, y)
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
        bfgs_in.BFGSop.solve(Hy, e)
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
        bfgs.parameters['H0inv'] = 'default'

        for ii in range(10):
            s.set_local((ii+1.0)*np.random.randn(dim))
            y.set_local(np.random.randn(dim))
            bfgs.BFGSop.update(s,y)

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
        bfgs.parameters['H0inv'] = 'default'
        bfgs.BFGSop.updated0 = False

        Hm = np.eye(dim)

        for ii in range(30):
            s.set_local((ii+1.0)*np.random.randn(dim))
            y.set_local(np.random.randn(dim))
            bfgs.BFGSop.update(s, y)

            # Update matrix form
            sk = s.array().reshape((dim,-1))
            yk = y.array().reshape((dim,-1))
            rho = bfgs.BFGSop.R[-1]
            Ir = np.eye(dim) - rho*sk.dot(yk.T)
            Hm = Ir.dot(Hm.dot(Ir.T)) + rho*sk.dot(sk.T)

        Hk = assemble_Hk(bfgs)

        diffHkHm = np.linalg.norm(Hk-Hm)/np.linalg.norm(Hk)
        print 'Diff between Hk and Hm = {}'.format(diffHkHm),
        if diffHkHm < 1e-12:   print '\t==>> OK!!'
        else:   print '\t\t\t<<<<==== PROBLEM!?'



def test4():
    """
    Check how quickly BFGS Hessian converges toward inverse Hessian of quadratic form
    f(x) = f(x0) + g^T (x-x0) + 1/2 (x-x0)^T B (x-x0)
    """
    def grad(x, x0, g, B):
        return g + B.dot(x-x0)

    mesh = dl.UnitSquareMesh(10, 10)
    V = dl.FunctionSpace(mesh, 'CG', 1)
    zp = ZeroPrior(V)
    model = Model([], zp, [])

    s, y, = dl.Vector(), dl.Vector()
    zp.init_vector(s, 0)
    zp.init_vector(y, 0)
    dim = len(s.array())
    print 'dim={}'.format(dim)

    for ii in range(7):
        print '\nTest {}'.format(ii)
        if ii == 0:
            # exact solution: Hk=I for all k provided H0=I
            B = np.eye(dim)
            invB = np.eye(dim)
            print 'B=Id'
        elif ii == 1:
            q, r = np.linalg.qr(np.random.randn(dim*dim).reshape((dim,dim)))
            l = np.linspace(1.,1.1,dim)
            print 'min(eig)={}, max(eig)={}'.format(np.min(l), np.max(l))
            B = q.dot(np.diag(l).dot(q.T))
            invB = q.dot(np.diag(1./l).dot(q.T))
        elif ii == 2:
            q, r = np.linalg.qr(np.random.randn(dim*dim).reshape((dim,dim)))
            l = np.linspace(1.,10.,dim)
            print 'min(eig)={}, max(eig)={}'.format(np.min(l), np.max(l))
            B = q.dot(np.diag(l).dot(q.T))
            invB = q.dot(np.diag(1./l).dot(q.T))
        elif ii == 3:
            q, r = np.linalg.qr(np.random.randn(dim*dim).reshape((dim,dim)))
            l = np.linspace(1.,1000.,dim)
            print 'min(eig)={}, max(eig)={}'.format(np.min(l), np.max(l))
            B = q.dot(np.diag(l).dot(q.T))
            invB = q.dot(np.diag(1./l).dot(q.T))
        elif ii == 4:
            q, r = np.linalg.qr(np.random.randn(dim*dim).reshape((dim,dim)))
            l = np.linspace(1e-1, 1.0, dim)
            print 'min(eig)={}, max(eig)={}'.format(np.min(l), np.max(l))
            B = q.dot(np.diag(l).dot(q.T))
            invB = q.dot(np.diag(1./l).dot(q.T))
        elif ii == 5:
            q, r = np.linalg.qr(np.random.randn(dim*dim).reshape((dim,dim)))
            l = np.linspace(1e-5, 1.0, dim)
            print 'min(eig)={}, max(eig)={}'.format(np.min(l), np.max(l))
            B = q.dot(np.diag(l).dot(q.T))
            invB = q.dot(np.diag(1./l).dot(q.T))
        else:
            q, r = np.linalg.qr(np.random.randn(dim*dim).reshape((dim,dim)))
            l = np.random.randn(dim)**2
            print 'min(eig)={}, max(eig)={}'.format(np.min(l), np.max(l))
            B = q.dot(np.diag(l).dot(q.T))
            invB = q.dot(np.diag(1./l).dot(q.T))
        print '|B*invB-I|={:.2e}, |invB*B-I|={:.2e}'.format(
        np.linalg.norm(B.dot(invB)-np.eye(dim)),
        np.linalg.norm(invB.dot(B)-np.eye(dim)))
        g = np.random.randn(dim).reshape((dim,-1))
        x0 = np.random.randn(dim).reshape((dim,-1))
        x1 = np.random.randn(dim).reshape((dim,-1))
        x2 = np.random.randn(dim).reshape((dim,-1))

        bfgs = BFGS(model)
        bfgs.parameters['H0inv'] = 'default'

        norminvB = np.linalg.norm(invB)

        for ii in range(1000):
            s.set_local(x2-x1)
            y.set_local(grad(x2,x0,g,B) - grad(x1,x0,g,B))
            bfgs.BFGSop.update(s, y)

            if ii%100 == 0:
                Hk = assemble_Hk(bfgs)
                diffBHk = np.linalg.norm(Hk-invB)/norminvB
                print 'ii={}; |Hk|={:.2e}, |inv(B)|={:.2e}, |Hk-inv(B)|/|inv(B)|={:.2e}'.format(
                ii, np.linalg.norm(Hk), norminvB, diffBHk)

            x1 = x2
            x2 = x1 + np.random.randn(dim).reshape((dim,-1))

        Hk = assemble_Hk(bfgs)
        diffBHk = np.linalg.norm(Hk-invB)/norminvB
        print 'ii={}; |Hk|={:.2e}, |inv(B)|={:.2e}, |Hk-inv(B)|/|inv(B)|={:.2e}'.format(
        ii, np.linalg.norm(Hk), norminvB, diffBHk)


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
    elif case == 4:
        test4()
    else:
        print 'Test case {} not implemented (yet)'.format(case)
        sys.exit(1)
