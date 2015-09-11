from dolfin import compile_extension_module, Vector, PETScKrylovSolver
import os
import numpy as np

def amg_method():
    S = PETScKrylovSolver()
    for pp in S.preconditioners():
        if pp[0] == 'ml_amg':
            return 'ml_amg'
        
    return 'petsc_amg'

abspath = os.path.dirname( os.path.abspath(__file__) )
sdir = os.path.join(abspath,"cpp_linalg")
header_file = open(os.path.join(sdir,"linalg.h"), "r")
code = header_file.read()
header_file.close()
cpp_sources = ["linalg.cpp"]  
cpp_module = compile_extension_module(
code=code, source_directory=sdir, sources=cpp_sources,
include_dirs=[".",  sdir])

def MatMatMult(A,B):
    s = cpp_module.cpp_linalg()
    return s.MatMatMult(A,B)

def MatPtAP(A,P):
    s = cpp_module.cpp_linalg()
    return s.MatPtAP(A,P)

def MatAtB(A,B):
    s = cpp_module.cpp_linalg()
    return s.MatAtB(A,B)

def to_dense(A):
    n  = A.size(0)
    m  = A.size(1)
    print n,m
    B = np.zeros( (n,m), dtype=np.float64)
    for i in range(0,n):
        [j, val] = A.getrow(i)
        B[i,j] = val
        
    return B

def getColoring(A,k):
    return cpp_module.Coloring(A,k)


def trace(A):
    n  = A.size(0)
    tr = 0.
    for i in range(0,n):
        [j, val] = A.getrow(i)
        tr += val[j == i]
    return tr

def get_diagonal(A, d, solve_mode=True):
    
    ej, xj = Vector(), Vector()

    if hasattr(A, "init_vector"):
        A.init_vector(ej,1)
        A.init_vector(xj,0)
    else:       
        A.get_operator().init_vector(ej,1)
        A.get_operator().init_vector(xj,0)
        
    ncol = ej.size()
    da = np.zeros(ncol, dtype=ej.array().dtype)
    
    for j in range(ncol):
        ej[j] = 1.
        if solve_mode:
            A.solve(xj, ej)
        else:
            A.mult(ej,xj)
        da[j] = xj[j]
        ej[j] = 0.
        
    d.set_local(da)

#def probe_diagonal_inv(Asparsity, Asolver, d, tol):
#    x, b, Ab = Vector(), Vector(), Vector()
#    
#    if hasattr(Asolver, "init_vector"):
#        Asolver.init_vector(x,1)
#        Asolver.init_vector(b,0)
#        Asolver.init_vector(Ab,0)
#    else:       
#        Asolver.get_operator().init_vector(x,1)
#        Asolver.get_operator().init_vector(b,0)
#        Asolver.get_operator().init_vector(Ab,0)
#        
#    b[0] = 1.
#    Asolver.mult(x,b)
#    err = 1.
#    
#    while err > tol:
#        Asparsity.mult(b,Ab)
        
        
    
    
def estimate_diagonal_inv_coloring(Asolver, coloring, d):
    x, b = Vector(), Vector()
    
    if hasattr(Asolver, "init_vector"):
        Asolver.init_vector(x,1)
        Asolver.init_vector(b,0)
    else:       
        Asolver.get_operator().init_vector(x,1)
        Asolver.get_operator().init_vector(b,0)
    
    ncolors = coloring.numberOfColors()
    print ncolors
    
    d.zero()
    for color in range(ncolors):
        coloring.markcolor(color, b, 1.)
        x.zero()
        Asolver.solve(x,b)
        coloring.x_dot_mult_b(color, x, b, d)


def estimate_diagonal_inv2(Asolver, k, d):
    x, b = Vector(), Vector()
    
    if hasattr(Asolver, "init_vector"):
        Asolver.init_vector(x,1)
        Asolver.init_vector(b,0)
    else:       
        Asolver.get_operator().init_vector(x,1)
        Asolver.get_operator().init_vector(b,0)
    
    num = np.zeros(b.array().shape, dtype = b.array().dtype)
    den = np.zeros(num.shape, dtype = num.dtype)
    for i in range(k):
        x.zero()
        b.set_local(np.random.randn(num.shape[0]))
        Asolver.solve(x,b)
        num = num +  ( x.array() * b.array() )
        den = den +  ( b.array() * b.array() )
        
    d.set_local( num / den )
        
def randn_perturb(x, std_dev):
    """Add a Gaussian random perturbation to x"""
    n = x.array().shape[0]
    noise = np.random.normal(0, 1, n)
    x.set_local(x.array() + std_dev*noise)