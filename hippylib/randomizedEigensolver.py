# Copyright (c) 2016, The University of Texas at Austin & University of
# California, Merced.
#
# All Rights reserved.
# See file COPYRIGHT for details.
#
# This file is part of the hIPPYlib library. For more information and source code
# availability see https://hippylib.github.io.
#
# hIPPYlib is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License (as published by the Free
# Software Foundation) version 3.0 dated June 2007.

from dolfin import Vector, MPI
from linalg import MultiVector, MatMvMult, MvDSmatMult, Solver2Operator
import numpy as np
import math
import matplotlib.pyplot as plt

"""
Randomized algorithms for the solution of Hermitian Eigenvalues Problems (HEP)
and Generalized Hermitian Eigenvalues Problems (GHEP).

In particular we provide an implementation of the single and double pass algorithms
and some convergence test.

REFERENCES:

Nathan Halko, Per Gunnar Martinsson, and Joel A. Tropp,
Finding structure with randomness:
Probabilistic algorithms for constructing approximate matrix decompositions,
SIAM Review, 53 (2011), pp. 217-288.

Arvind K. Saibaba, Jonghyun Lee, Peter K. Kitanidis,
Randomized algorithms for Generalized Hermitian Eigenvalue Problems with application
to computing Karhunen-Loeve expansion,
Numerical Linear Algebra with Applications, to appear.

NOTE: This routines are only serial!!
"""

def singlePass(A,Omega,k):
    """
    The single pass algorithm for the HEP as presented in [1].
    Inputs:
    - A: the operator for which we need to estimate the dominant eigenpairs.
    - Omega: a random gassian matrix with m >= k columns.
    - k: the number of eigenpairs to extract.
    
    Outputs:
    - d: the estimate of the k dominant eigenvalues of A
    - U: the estimate of the k dominant eigenvectors of A. U^T U = I_k
    """
    w = Vector()
    y = Vector()
    A.init_vector(w,1)
    A.init_vector(y,0)
    
    mpi_comm = w.mpi_comm()
    nprocs = MPI.size(mpi_comm)
    if nprocs > 1:
        raise Exception("function singlePass is only serial")
    
    nvec  = Omega.shape[1]
    
    assert(nvec >= k )
    
    Y = np.zeros(Omega.shape)
    
    for ivect in range(0,nvec):
        w.set_local(Omega[:,ivect])
        A.mult(w,y)
        Y[:,ivect] = y.array()
                
    Q,_ = np.linalg.qr(Y)
        
    Zt = np.dot(Omega.T, Q)
    Wt = np.dot(Y.T, Q)
        
    Tt = np.linalg.solve(Zt, Wt)
                
    T = .5*Tt + .5*Tt.T
        
    d, V = np.linalg.eigh(T)
    sort_perm = d.argsort()
        
    sort_perm = sort_perm[::-1]
    d = d[sort_perm[0:k]]
    V = V[:, sort_perm[0:k]] 
        
    U = np.dot(Q,V)
        
    return d, U

def doublePass(A,Omega,k):
    """
    The double pass algorithm for the HEP as presented in [1].
    Inputs:
    - A: the operator for which we need to estimate the dominant eigenpairs.
    - Omega: a random gassian matrix with m >= k columns.
    - k: the number of eigenpairs to extract.
    
    Outputs:
    - d: the estimate of the k dominant eigenvalues of A
    - U: the estimate of the k dominant eigenvectors of A. U^T U = I_k
    """
    w = Vector()
    y = Vector()
    A.init_vector(w,1)
    A.init_vector(y,0)
    
    mpi_comm = w.mpi_comm()
    nprocs = MPI.size(mpi_comm)
    if nprocs > 1:
        raise Exception("function doublePass is only serial")
    
    nvec  = Omega.shape[1]
    
    assert(nvec >= k )
    
    Y = np.zeros(Omega.shape)
    
    for ivect in range(0,nvec):
        w.set_local(Omega[:,ivect])
        A.mult(w,y)
        Y[:,ivect] = y.array()
                
    Q,_ = np.linalg.qr(Y)
    
    AQ = np.zeros(Omega.shape)
    for ivect in range(0,nvec):
        w.set_local(Q[:,ivect])
        A.mult(w,y)
        AQ[:,ivect] = y.array()
                
    T = np.dot(Q.T, AQ)
        
    d, V = np.linalg.eigh(T)
    sort_perm = d.argsort()
        
    sort_perm = sort_perm[::-1]
    d = d[sort_perm[0:k]]
    V = V[:, sort_perm[0:k]] 
        
    U = np.dot(Q,V)
        
    return d, U

def singlePassG(A, B, Binv, Omega,k, check_Bortho = False, check_Aortho=False, check_residual = False):
    """
    The single pass algorithm for the GHEP as presented in [2].
    B-orthogonalization is achieved using the PreCholQR algorithm.
    
    Inputs:
    - A: the operator for which we need to estimate the dominant generalized eigenpairs.
    - B: the rhs operator
    - Omega: a random gassian matrix with m >= k columns.
    - k: the number of eigenpairs to extract.
    
    Outputs:
    - d: the estimate of the k dominant eigenvalues of A
    - U: the estimate of the k dominant eigenvectors of A. U^T B U = I_k
    """
    w = Vector()
    ybar = Vector()
    y = Vector()
    A.init_vector(w,1)
    A.init_vector(y,0)
    A.init_vector(ybar,0)
    
    mpi_comm = w.mpi_comm()
    nprocs = MPI.size(mpi_comm)
    if nprocs > 1:
        raise Exception("function singlePassG is only serial")
    
    nvec  = Omega.shape[1]
    
    assert(nvec >= k )
    
    Ybar = np.zeros(Omega.shape)
    Y = np.zeros(Omega.shape)
    
    for ivect in range(0,nvec):
        w.set_local(Omega[:,ivect])
        A.mult(w,ybar)
        Binv.solve(y, ybar)
        Ybar[:,ivect] = ybar.array()
        Y[:,ivect] = y.array()
                
    Z,_ = np.linalg.qr(Y)
    BZ = np.zeros(Omega.shape)
    for ivect in range(0,nvec):
        w.set_local(Z[:,ivect])
        B.mult(w,y)
        BZ[:, ivect] = y
        
    R = np.linalg.cholesky( np.dot(Z.T,BZ ))
    Q = np.linalg.solve(R, Z.T).T
    BQ = np.linalg.solve(R, BZ.T).T
    
    Xt = np.dot(Omega.T, BQ)
    Wt = np.dot(Ybar.T, Q)
    Tt = np.linalg.solve(Xt,Wt)
                
    T = .5*Tt + .5*Tt.T
        
    d, V = np.linalg.eigh(T)
    sort_perm = d.argsort()
        
    sort_perm = sort_perm[::-1]
    d = d[sort_perm[0:k]]
    V = V[:, sort_perm[0:k]] 
        
    U = np.dot(Q,V)
    
    if check_Bortho:
        BorthogonalityTest(B, U)
    if check_Aortho: 
        AorthogonalityCheck(A, U, d)
    if check_residual:
        residualCheck(A,B, U, d)
        
    return d, U

def doublePassG(A, B, Binv, Omega, k, s = 1, check = False):
    """
    The double pass algorithm for the GHEP as presented in [2].
    
    Inputs:
    - A: the operator for which we need to estimate the dominant generalized eigenpairs.
    - B: the rhs operator
    - Omega: a random gassian matrix with m >= k columns.
    - k: the number of eigenpairs to extract.
    
    Outputs:
    - d: the estimate of the k dominant eigenvalues of A
    - U: the estimate of the k dominant eigenvectors of A. U^T B U = I_k
    """        
    nvec  = Omega.nvec()
    
    assert(nvec >= k )
    
    Ybar = MultiVector(Omega[0], nvec)
    Q = MultiVector(Omega)
    for i in range(s):
        MatMvMult(A, Q, Ybar)
        MatMvMult(Solver2Operator(Binv), Ybar, Q)
    
    Q.Borthogonalize(B)
    AQ = MultiVector(Omega[0], nvec)
    MatMvMult(A, Q, AQ)
    
    T = AQ.dot_mv(Q)
                        
    d, V = np.linalg.eigh(T)
    sort_perm = d.argsort()
        
    sort_perm = sort_perm[::-1]
    d = d[sort_perm[0:k]]
    V = V[:, sort_perm[0:k]] 
        
    U = MultiVector(Omega[0], k)
    MvDSmatMult(Q, V, U)
    
    if check:
        checkAll(A,B, U, d)
            
    return d, U


def checkAll(A,B, U, d):
    nvec  = U.nvec()
    AU = MultiVector(U[0], nvec)
    BU = MultiVector(U[0], nvec)
    MatMvMult(A, U, AU)
    MatMvMult(B, U, BU)
    
    # Residual checks
    diff = MultiVector(AU)
    diff.axpy(-d, BU)
    res_norms = diff.norm("l2")
    
    # B-ortho check
    UtBU = BU.dot_mv(U)
    err = UtBU - np.eye(nvec, dtype=UtBU.dtype)
    err_Bortho = np.linalg.norm(err, 'fro')
    
    #A-ortho check
    V = MultiVector(U)
    scaling = np.power(d, -0.5)
    V.scale(scaling)
    AU.scale(scaling)
    VtAV = AU.dot_mv(V)
    err = VtAV - np.eye(nvec, dtype=VtAV.dtype)
    err_Aortho = np.linalg.norm(err, 'fro')
    
    mpi_comm = U[0].mpi_comm()
    rank = MPI.rank(mpi_comm)
    if rank == 0:
        print "|| UtBU - I ||_F = ", err_Bortho
        print "|| VtAV - I ||_F = ", err_Aortho, " with V = U D^{-1/2}"
        print "lambda", "||Au - lambdaBu||_2"
        for i in range(res_norms.shape[0]):
            print "{0:5e} {1:5e}".format(d[i], res_norms[i])
    
    
def BorthogonalityTest(B, U):
    """
    Test the frobenious norm of  U^TBU - I_k
    """
    BU = MultiVector(U[0], U.nvec())
    nvec  = U.nvec()
    
    MatMvMult(B, U, BU)        
    UtBU = BU.dot_mv(U)
    err = UtBU - np.eye(nvec, dtype=UtBU.dtype)
    err_norm = np.linalg.norm(err, 'fro')
    
    mpi_comm = U[0].mpi_comm()
    rank = MPI.size(mpi_comm)
    if rank == 0:
        print "|| UtBU - I ||_F = ", err_norm
    
def AorthogonalityCheck(A, U, d):
    """
    Test the frobenious norm of  D^{-1}(U^TAU) - I_k
    """
    nvec  = U.nvec()
    V = MultiVector(U)
    scaling = np.power(d, -0.5)
    V.scale(scaling)
    AV = MultiVector(U[0], U.nvec())
    MatMvMult(A, V, AV)
    
    VtAV = AV.dot_mv(V)
    err = VtAV - np.eye(nvec, dtype=VtAV.dtype)
    err_norm = np.linalg.norm(err, 'fro')

    
    mpi_comm = U[0].mpi_comm()
    rank = MPI.rank(mpi_comm)
    if rank == 0:
        print "|| VtAV - I ||_F = ", err_norm, " with V = U D^{-1/2}"
            
def residualCheck(A,B, U, d):
    """
    Test the l2 norm of the residual:
    r[:,i] = d[i] B U[:,i] - A U[:,i]
    """
    AU = MultiVector(U[0], U.nvec())
    BU = MultiVector(U[0], U.nvec())
    MatMvMult(A, U, AU)
    MatMvMult(B, U, BU)
    AU.axpy(-d, BU)
    norms = AU.norm("l2")
    
    mpi_comm = U[0].mpi_comm()
    rank = MPI.rank(mpi_comm)
    
    if rank == 0:    
        print "lambda", "||Au - lambdaBu||"
        for i in range(norms.shape[0]):
            print d[i], norms[i]
        