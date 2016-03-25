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

## Compute direct and adjoint sensitivities for 
## elliptic partial differential equation
##
##   min  1/2 * ||u - ud||^2
##    m
##   where u is the solution of
##
##   - div (m * grad u) = f     on Omega
##                    u = 0     on bdry(Omega)
##
## for given force f, and data ud.
## The data ud is constructed from u by adding noise.
## and consider m = sum_{i=1}^7 m_i phi_i(x), where
## phi_1(x) = 1; phi_2(x) = sin(2*pi*x); phi_3(x) = sin(2*pi*y), etc.

# Import dependencies
from dolfin import *
import numpy as np

# create mesh and define function spaces
nx = 32
ny = 32
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, 'Lagrange', 1)

# define function for state and adjoint
u = Function(V)
p = Function(V)

# define Trial and Test Functions
du_trial, p_trial, u_trial = TrialFunction(V), TrialFunction(V), TrialFunction(V)
du_test,  p_test,  u_test  = TestFunction(V), TestFunction(V), TestFunction(V)

# initialize input functions
f  = Constant("1.0")
u0 = Constant("0.0")
m  = Constant("2.0")

# set up dirichlet boundary conditions
def u0_boundary(x,on_boundary):
    return on_boundary
bc = DirichletBC(V, u0, u0_boundary)

# variational forms
var_state     = inner( m * nabla_grad(u_trial), nabla_grad(u_test)) * dx
var_rhs_state = f * u_test * dx
var_sens      = inner( m * nabla_grad(du_trial), nabla_grad(du_test)) * dx

# assemble just for getting the size
mat_state, rhs_state = assemble_system(var_state, var_rhs_state, bc)

# solve state equation
u = Function(V,name="state")
solve(var_state == var_rhs_state, u, bc)
File("state.pvd") << u

## set up misfit
ud = Function(V,name="data")
ud.assign(u)
## perturb state solution and create synthetic measurements ud
## ud = u + ||u||/SNR * random.normal
noise_level = 0.05
MAX         = ud.vector().norm("linf")
noise       = Vector()
mat_state.init_vector(noise,1)
noise.set_local( noise_level * MAX * np.random.normal(0, 1, len(ud.vector().array())) )
#bc.apply(noise)
ud.vector().axpy(1., noise)
File("data.pvd") << ud

var_misfit = inner( (u - ud), u_test)*dx
vec_misfit = assemble(var_misfit)

# solve adjoint problem
var_adj = inner( m * nabla_grad(p_trial), nabla_grad(p_test)) * dx
rhs_adj = -var_misfit
p = Function(V,name="adj")
solve(var_adj==rhs_adj,p,bc)
File("adj.pvd") << p

##################################
# first component of the gradient
##################################
dm = Expression('1.0')
var_rhs_sens  = inner( -dm * nabla_grad(u), nabla_grad(du_test) ) * dx
du = Function(V,name="sens1")
solve(var_sens==var_rhs_sens,du,bc)
File("sens1.pvd") << du
grad_sens1 = vec_misfit.inner(du.vector())
print grad_sens1

# compute the sensitivity using the assemble (solve mat_sense du = rhs_sens)
var_rhs_sens  = inner( -dm * nabla_grad(u), nabla_grad(du_test) ) * dx
mat_sens, rhs_sens = assemble_system(var_sens, var_rhs_sens, bc)
solve(mat_sens, du.vector(), rhs_sens)
grad_fd_sens1 = vec_misfit.inner(du.vector())
print grad_fd_sens1

# adjoint approach
grad_adj1 = assemble(inner (dm * nabla_grad(u), nabla_grad(p) )* dx)
print grad_adj1

##################################
# second component of the gradient
##################################
# sensitivity approach
dm = Expression("sin(2*pi*x[1])")
var_rhs_sens  = inner( -dm * nabla_grad(u), nabla_grad(du_test) ) * dx
du = Function(V,name="sens2")
solve(var_sens==var_rhs_sens,du,bc)
File("sens2.pvd") << du
grad_sens2 = vec_misfit.inner(du.vector())
print grad_sens2

# adjoint approach
grad_adj2 = assemble(inner (dm * nabla_grad(u), nabla_grad(p) )* dx)
print grad_adj2
