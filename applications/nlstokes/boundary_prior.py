'''
Created on Oct 8, 2015

@author: uvilla
'''
from dolfin import *
import sys
from dolfin.cpp.la import PETScKrylovSolver
sys.path.append( "../../" )
import numpy as np
from hippylib import *

nx = 16
ny = 16
nz = 16

mesh = Mesh("sphere.xml")
Vh = FunctionSpace(mesh, "CG", 1)

uh = TrialFunction(Vh)
vh = TestFunction(Vh)

class MyBoundary(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary
    
boundary_markers = FacetFunction("size_t",mesh)
boundary_markers.set_all(0)
Gamma = MyBoundary()
Gamma.mark(boundary_markers,1)
dss = Measure("ds")[boundary_markers]

n = FacetNormal(mesh)

fExp = Expression("10*exp( -100*(x[0]-1)*(x[0]-1) - 100*x[1]*x[1] - 100*x[2]*x[2])")
avarf = uh*vh*dss(1) + 0.1*inner(grad(uh) - dot(outer(n,n),grad(uh)), grad(vh) - dot(outer(n,n),grad(vh)) )*dss(1)
#avarf = uh*vh*dss(1) + 0.1*inner(grad(uh), grad(vh) )*dss(1)
fvarf = fExp*vh*dss(1)

A = assemble(avarf,keep_diagonal=True)
A.ident_zeros()

rhs = assemble(fvarf)

solver= PETScKrylovSolver("cg", "amg")
solver.set_operator(A)

sol = Function(Vh)

solver.solve(sol.vector(), rhs)

plot(sol, interactive=True)