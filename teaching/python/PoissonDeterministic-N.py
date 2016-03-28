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


## Newton solution of a coefficient field inversion in an
## elliptic partial differential equation
##
##   min  1/2 * ||u - uobs||^2 + gamma/2 * ||grad a||^2,
##    a
##        where u is the solution of
##
##             - div (a * grad u) = f     on Omega
##                              u = 0     on bdry(Omega)
##
## for given force f, gamma >= 0 and data uobs.
## The data uobs is constructed using a "true" parameter field atrue

from dolfin import *

import numpy as np
import matplotlib.pyplot as plt
import time
import logging

start = time.clock()

logging.getLogger('FFC').setLevel(logging.WARNING)

set_log_active(False)
print "Messages and Warnings Disabled"
set_log_level(DEBUG)
np.random.seed(seed=1)

def Norm_Check(*args):
    for arg in args:
        var = eval(arg)
        if type(var) == Function:
            norm_var = np.linalg.norm(var.vector().array())
        elif type(var) == Matrix or type(var) == GenericVector:
            norm_var = np.linalg.norm(var.array())
        else:
            norm_var = 0
            print "Norm %s:" % (arg), norm_var

# define cost function
def cost(u, ud, a, W, R):
    diff = u.vector() - ud.vector()
    reg = 0.5 * a.vector().inner(R*a.vector() ) 
    misfit = 0.5 * diff.inner(W * diff)
    return [reg + misfit, misfit, reg]

# define (Gauss-Newton) Hessian apply H * v
def Hess_GN (v, R, C, A, W):
    solve (A, du, - (C * v))
    solve (A, dp, - (W * du))
    CT_dp = Vector()
    C.init_vector(CT_dp, 1)
    C.transpmult(dp, CT_dp)
    H_V = R * v + CT_dp
    return H_V

# define (Newton) Hessian apply H * v
def Hess_Newton (v, R, C, A, W, Wua):
    RHS = -(C * v)
    bc2.apply(RHS)
    solve (A, du, RHS)
    RHS = -(W * du) -  Wua * v
    bc2.apply(RHS)
    solve (A, dp, RHS)
    CT_dp = Vector()
    C.init_vector(CT_dp, 1)
    C.transpmult(dp, CT_dp)
    Wua_du = Vector()
    Wua.init_vector(Wua_du, 1)
    Wua.transpmult(du, Wua_du)
    H_V = R*v + CT_dp + Wua_du
    return H_V

# Creat Class MyLinearOperator to perform Hessian function
class MyLinearOperator(LinearOperator):
    cgiter = 0
    def __init__(self, R, C, A, W, Wua):
        LinearOperator.__init__(self, a_delta, a_delta) # why this line?
        self.R = R
        self.C = C
        self.A = A
        self.W = W
        self.Wua = Wua

    # Hessian performed on x, output as generic vector y
    def mult(self, x, y):
        self.cgiter += 1
        if iter <= 500:
            y.set_local(Hess_GN (x, self.R, self.C, self.A, self.W).array() )
        else:
            y.set_local(Hess_Newton (x, self.R, self.C, self.A, self.W, self.Wua).array() )

# define parameters
tol = 1e-8
gamma = 5e-11
c = 1e-4
maxiter = 100
eps = 1e-4

# initialize iter counters
iter = 0
total_cg_iter = 0
solution = 0

# create mesh and define function spaces
nx = 64
ny = 64
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, 'Lagrange', 1)
V2 = FunctionSpace(mesh, 'Lagrange', 2)

# define Trial and Test Functions
u, p, ud, g, a_delta, a_try     = TrialFunction(V2), TrialFunction(V2), TrialFunction(V2), TrialFunction(V), TrialFunction(V), TrialFunction(V)
u_test, p_test, ud_test, g_test = TestFunction(V2), TestFunction(V2), TestFunction(V2), TestFunction(V)
p = Function(V2)

# initialize input functions
atrue = interpolate(Expression('8. - 4.*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2) < pow(0.2,2))'), V)
f = Constant("1.0")
u0 = Constant("0.0")
a = interpolate(Expression("4."),V)

# set up dirichlet boundary conditions
def u0_boundary(x,on_boundary):
    return on_boundary
bc = DirichletBC(V, u0, u0_boundary)
bc2 = DirichletBC(V2, u0, u0_boundary)

# weak form for setting up the synthetic observations
a_goal = inner( atrue * nabla_grad(ud), nabla_grad(ud_test)) * dx
L_goal = f * ud_test * dx

# solve the forward/state problem
goal_A, goal_b = assemble_system(a_goal, L_goal, bc2)
ud = Function(V2)
solve(goal_A, ud.vector(), goal_b)

# perturb state solution and create synthetic measurements ud
MAX = ud.vector().norm("linf")
noise = Vector()
goal_A.init_vector(noise,1)
noise.set_local( .01 * MAX * np.random.normal(0, 1, len(ud.vector().array())) )
ud.vector().axpy(1., noise)
bc2.apply(ud.vector())

# weak form for setting up the state equation
a_state = inner( a * nabla_grad(u), nabla_grad(u_test)) * dx
L_state = f * u_test * dx
W_equ   = inner(u, u_test) * dx

# # weak form for setting up the right hand side of the adjoint
u = Function(V2)
L_adjoint = -inner(u - ud, u_test) * dx

# weak form for setting up matrices
Wua_equ = inner(a_delta * nabla_grad(p_test), nabla_grad(p)) * dx
C_equ   = inner( a_delta * nabla_grad(u), nabla_grad(u_test)) * dx
M_equ   = inner(g, g_test) * dx
R_equ   = gamma * inner(nabla_grad(g), nabla_grad(g_test)) * dx

# assemble matrices M, W, and R
M = assemble(M_equ)
W = assemble(W_equ)
R = assemble(R_equ)

# solve state equation
A, state_b = assemble_system (a_state, L_state, bc2)
solve (A, u.vector(), state_b)

# initializations
g, a_delta = Vector(), Vector()
R.init_vector(a_delta,0)
R.init_vector(g,0)

du, dp = Vector(), Vector()
W.init_vector(du,1)
W.init_vector(dp,0)

a_prev, a_diff = Function(V), Function(V)

print "Nit  CGit  cost          misfit        reg           sqrt(-G*D)    ||grad||       alpha  tolcg"

while iter <  maxiter and solution == 0:
    iter += 1

    # assemble matrix C
    C =  assemble(C_equ)

    # solve the adoint problem
    A, adjoint_RHS = assemble_system(a_state, L_adjoint, bc2)
    solve(A, p.vector(), adjoint_RHS)

    # assemble W and R
    Wua = assemble (Wua_equ)

    # evaluate the  gradient
    CT_p = Vector()
    C.init_vector(CT_p,1)
    C.transpmult(p.vector(), CT_p)
    MG = CT_p + R * a.vector()
    solve(M, g, MG)

    # calculate the norm of the gradient
    grad2 = g.inner(MG)
    gradnorm = sqrt(grad2)

    # set the CG tolerance
    if iter == 1:
        gradnorm_ini = gradnorm
    tolcg = min(0.5, (gradnorm/gradnorm_ini))

    # define the Hessian apply operator (with preconditioner)
    Hess_Apply = MyLinearOperator(R, C, A, W, Wua )
    P = R + 100 * gamma * M
    solver = PETScKrylovSolver("cg", "amg")
    solver.set_operators(Hess_Apply, P)
    solver.parameters["relative_tolerance"] = tolcg

    # solve the Newton system H a_delta = - MG
    solver.solve(a_delta, -MG)
    total_cg_iter += Hess_Apply.cgiter

    if iter == 1:   # first iteration cost (before update of a)
        [cost_old, misfit_old, reg_old] = cost(u, ud, a, W, R)

    # linesearch
    alpha = 1
    descent = 0
    no_backtrack = 0
    a_prev.assign(a)
    while descent == 0 and no_backtrack < 10:
        a.vector().axpy(alpha, a_delta )

        # solve the state/forward problem
        state_A, state_b = assemble_system(a_state, L_state, bc2)
        solve(state_A, u.vector(), state_b)

        # evaluate cost
        [cost_new, misfit_new, reg_new] = cost(u, ud, a, W, R)

        # check if Armijo conditions are satisfied
        if cost_new < cost_old - alpha * c * grad2:
            cost_old = cost_new
            descent = 1
        else:
            no_backtrack += 1
            alpha *= 0.5
            a.assign(a_prev)  # reset a

    # calculate ||a-atrue|| and sqrt(-G * D)
    a_diff.assign(a)
    a_diff.vector().axpy(-1., interpolate(atrue,V).vector() )
    graddir = sqrt(- MG.inner(a_delta) )

    sp = ""
    print "%1d %2s %2d %3s %8.5e %1s %8.5e %1s %8.5e %1s %8.5e %1s %8.5e %1s %5.2f %1s %5.3e" % \
        (iter, sp, Hess_Apply.cgiter, sp, cost_new, sp, misfit_new, sp, reg_new, sp, \
         graddir, sp, gradnorm, sp, alpha, sp, tolcg)

    # check for convergence
    if sqrt(grad2) < tol and iter > 1:
        solution = 1
        print "Newton's method converged in ",iter,"  iterations"
        print "Total number of CG iterations: ", total_cg_iter
if solution == 0:
    print "Newton's method did not converge in ", maxiter, " iterations"

# Dump solution to file in VTK format
u.rename("u","ignore_this")
File('poisson_state.pvd') << u

p.rename("p","ignore_this")
File('poisson_adj.pvd') << p

a.rename("a","ignore_this")
File("poisson_control.pvd") << a

a.assign(interpolate(atrue,V))
a.rename("atrue","ignore_this")
File("poisson_control_true.pvd") << a
