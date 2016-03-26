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
# Software Foundation) version 2.1 dated February 1999.

import dolfin as dl
import os
import sys
sys.path.append( "../")
from hippylib import *

abspath = os.path.dirname( os.path.abspath(__file__) )
sdir = os.path.join(abspath,"cpp_coloring")
header_file = open(os.path.join(sdir,"coloring.h"), "r")
code = header_file.read()
header_file.close()
cpp_sources = ["coloring.cpp"]  
cpp_module = dl.compile_extension_module(
code=code, source_directory=sdir, sources=cpp_sources,
include_dirs=[".",  sdir])

def getColoring(A,k):
    """
    Apply a coloring algorithm the adjacency graph of A^k.
    """
    return cpp_module.Coloring(A,k)

def estimate_diagonal_inv_coloring(Asolver, coloring, d):
    """
    Use a probing algorithm to estimate the diagonal of A^{-1}.
    - Asolver:  a linear solver for the operator A.
    - coloring: a coloring based on the adjacency graph of A^k,
                for some power k. The number of linear system to
                be solved is equal to the number of colors.
    - d:        the estimated diagonal of A^{-1}.
    
    REFERENCE:
    Jok M Tang and Yousef Saad,
    A probing method for computing the diagonal of a matrix inverse,
    Numerical Linear Algebra with Applications, 19 (2012), pp. 485-501
    """
    
    x, b = dl.Vector(), dl.Vector()
    
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

def probingEstimator(prior, path_len):
    """
    We use a probing algorithm to approximate the diagonal entries of R^{-1}.
    The number of linear system in R to solve depends on the sparsity pattern of R.
    Even if the computational cost is much lower than the exact method, the complexity still increase linearly
    with the size of the problem. path_len represent the power of R used to compute the adjacency matrix for the
    graph coloring algorithm. See function estimate_diagonal_inv_coloring for details
    """
    pw_var = dl.Vector()
    prior.init_vector(pw_var,0)
    if type(prior.R) == dl.Matrix:
        coloring = getColoring(prior.R, path_len)
    else:
        Rsparsity = MatPtAP(prior.M, prior.A)
        coloring = getColoring(Rsparsity, path_len)
                
    estimate_diagonal_inv_coloring(prior.Rsolver, coloring, pw_var)
