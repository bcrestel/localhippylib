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

# Coefficient field inversion in an elliptic partial differential equation
# 
# Consider the following problem:
#
# min_a J(a):=1/2 int_Omega (u-ud)^2 dx +gamma/2 int_Omega | grad a|^2 dx
# 
# where u is the solution of
#
# -div (exp{a} grad u) + v grad u = f in Omega,
#               u = 0 on partial Omega.
# 
# Here a  the unknown coefficient field, ud denotes (possibly noisy) data, $f\in H^{-1}(Omega)$ a given force, and $ gamma >= 0$ the regularization parameter.

###############################################################################
# 0. Import dolfin and other modules
from dolfin import *

import math
import numpy as np
import time
import logging

import matplotlib.pyplot as plt
import nb

from cgsolver import CGSolver
from cgsolverSteihaug import CGSolverSteihaug

start = time.clock()

logging.getLogger('FFC').setLevel(logging.WARNING)
logging.getLogger('UFL').setLevel(logging.WARNING)
set_log_active(False)

np.random.seed(seed=1)

###############################################################################
# 1. Create the mesh and define the finite element spaces
nx = 80
ny = 80
mesh = UnitSquareMesh(nx, ny)
Va = FunctionSpace(mesh, 'Lagrange', 1)
Vu = FunctionSpace(mesh, 'Lagrange', 2)

###############################################################################
# 2. Define the true parameter and an initial guess for the parameter to be reconstructed
atrue = interpolate(Expression('log(2 + 6*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2) > pow(0.2,2)) )'),Va)
a = interpolate(Expression("log(8.0)"),Va)

plt.figure(figsize=(15,5))
nb.plot(mesh,subplot_loc=121, mytitle="Mesh", show_axis='on')
nb.plot(atrue,subplot_loc=122, mytitle="True parameter field")

###############################################################################
# 3. Define the finite element functions for state and adjoint variables
u = Function(Vu)
p = Function(Vu)

###############################################################################
# 4. Define the trial and test functions used by FEniCS
#    to assemble the weak form of the first and second variations
u_trial, p_trial, a_trial = TrialFunction(Vu), TrialFunction(Vu), TrialFunction(Va)
u_test, p_test, a_test = TestFunction(Vu), TestFunction(Vu), TestFunction(Va)

###############################################################################
# 5. Initialize the parameters of the PDE and the boundary conditions
f = Constant("10.0")
u0 = Constant("0.0")
v = Expression(("30", "0"))

def boundary(x,on_boundary):
    return on_boundary

bc_state = DirichletBC(Vu, u0, boundary)
bc_adj = DirichletBC(Vu, Constant(0.), boundary)


###############################################################################
# 6. Generate the synthetic observations by solving the forward problem.
#    Remember to modify the definition of a_goal when you solve the
#    advection diffusion problem.

# Set the noise level 
# When you change the noise level, remember to also change the regularization parameter beta
# according to the discrepancy criterion.
if 0:
    noise_level = 0.01 #0.02
    beta = 5e-8        #2e-8
else:
    noise_level = 0.1
    beta = 1.e-6        

# Define weak form of the forward problem for setting up the synthetic observations
a_goal = inner(exp(atrue) * nabla_grad(u_trial), nabla_grad(u_test)) * dx +\
         inner(v, nabla_grad(u_trial))*u_test * dx
L_goal = f * u_test * dx

# Solve the forward/state problem to generate synthetic observations
goal_A, goal_b = assemble_system(a_goal, L_goal, bc_state)

utrue = Function(Vu)
solve(goal_A, utrue.vector(), goal_b)

ud = Function(Vu)
ud.assign(utrue)

# Perturb state solution and create synthetic measurements ud
# ud = u + eta,
# where eta_i ~ N(0, sigma^2) and sigma = noise_level * || u_true ||_inf
MAX = ud.vector().norm("linf")
noise = Vector()
goal_A.init_vector(noise,1)
noise.set_local( noise_level * MAX * np.random.normal(0, 1, len(ud.vector().array())) )
bc_adj.apply(noise)

ud.vector().axpy(1., noise)

print "Using noise level = ", noise_level
print "(Noise variance)/2 = ", math.pow(noise_level*MAX, 2)/2

nb.multi1_plot([utrue, ud], ["State solution with atrue", "Synthetic observations"])

###############################################################################
# 7. Precompute the finite element matrices W and R for the discretization of
#    the cost functional.
#    Define the function cost that computes the total cost, the misfit componenent,
#    and the regularization component for a given state u and parameter a
#    

W_equ   = inner(u_trial, u_test) * dx
R_equ   = Constant(beta) * inner(nabla_grad(a_trial), nabla_grad(a_test)) * dx

W = assemble(W_equ)
R = assemble(R_equ)

def cost(u, ud, a, W, R):
    diff = u.vector() - ud.vector()
    reg = 0.5 * a.vector().inner(R*a.vector() ) 
    misfit = 0.5 * diff.inner(W * diff)
    return [reg + misfit, misfit, reg]

###############################################################################
# 8. Define the linear and bilinear form the the state/adjoint equations, and
#    incremental forward/adjoint equations.
#    You will need to modify the implementation of a_state and a_adj when you solve
#    the advection-diffusion inversion problem

# weak form for setting up the state equation
a_state = inner(exp(a) * nabla_grad(u_trial), nabla_grad(u_test)) * dx +\
          inner(v, nabla_grad(u_trial))*u_test * dx
L_state = f * u_test * dx

# weak form for setting up the adjoint equation
a_adj = inner(exp(a) * nabla_grad(p_trial), nabla_grad(p_test)) * dx -\
        inner(v, nabla_grad(p_trial))*p_test * dx
L_adj = -inner(u - ud, p_test) * dx

# additional weak form for the incremental fwd/adj
Wua_equ = inner(exp(a) * a_trial * nabla_grad(p_test), nabla_grad(p)) * dx
C_equ   = inner(exp(a) * a_trial * nabla_grad(u), nabla_grad(u_test)) * dx
Raa_equ = inner(exp(a) * a_trial * a_test *  nabla_grad(u),  nabla_grad(p)) * dx

M_equ   = inner(a_trial, a_test) * dx

# assemble matrix M
M = assemble(M_equ)

###############################################################################
# 9. Solve state equation for the initial guess of the parameter a and evaluate
#    the cost functional
state_A, state_b = assemble_system (a_state, L_state, bc_state)
solve (state_A, u.vector(), state_b)

[cost_old, misfit_old, reg_old] = cost(u, ud, a, W, R)

plt.figure(figsize=(15,5))
nb.plot(a,subplot_loc=121, mytitle="a_ini", vmin=atrue.vector().min(), vmax=atrue.vector().max())
nb.plot(u,subplot_loc=122, mytitle="u(a_ini)")

###############################################################################
# 10. Define the (Gauss-Newton) Hessian apply operator.
#     The method mult(v, Hv) will apply the (Gauss-Newton) Hessian to the vector v:
#     Hv = H*v
#     When using the conjugate gradient method, remember that while the Gauss-Newton Hessian
#     is always positive defined, the Hessian may be undefine far away from the minimum. 
#     
#     We also define a preconditioner P for the (Gauss-Newton) Hessian: P = R + beta M.

class GaussNewtonHessian:
    def __init__(self, R, C, A, adj_A, W):
        self.R = R
        self.C = C
        self.A = A
        self.adj_A = adj_A
        self.W = W
        
    def init_vector(self,x,dim):
        return R.init_vector(x,dim)
        
    def mult(self, v, Hv):
        rhs = -(self.C * v)
        bc_adj.apply(rhs)
        solve (self.A, du, rhs)
        rhs = - (self.W * du)
        bc_adj.apply(rhs)
        solve (self.adj_A, dp, rhs)
        self.C.transpmult(dp, Hv)
        Hv.axpy(1., self.R*v)
        
class Hessian(LinearOperator):
    def __init__(self, R, Raa, C, A, adj_A, W, Wua):
        LinearOperator.__init__(self, a_delta, a_delta)
        self.R = R
        self.Raa = Raa
        self.C = C
        self.A = A
        self.adj_A = adj_A
        self.W = W
        self.Wua = Wua
        
    def init_vector(self,x,dim):
        return R.init_vector(x,dim)
        
    def mult(self,v, Hv):
        rhs = -(self.C * v)
        bc_adj.apply(rhs)
        solve (self.A, du, rhs)
        rhs = -(self.W * du) -  self.Wua * v
        bc_adj.apply(rhs)
        solve (self.adj_A, dp, rhs)
        self.C.transpmult(dp, Hv)
        Wua_du = Vector()
        self.Wua.init_vector(Wua_du, 1)
        self.Wua.transpmult(du, Wua_du)
        Hv.axpy(1., Wua_du)
        Hv.axpy(1., self.R*v)
        Hv.axpy(1., self.Raa*v)
        
        
P = R + beta * M
prec = PETScKrylovSolver("cg", "amg")
prec.set_operator( P)
prec.parameters["relative_tolerance"] = 1e-9
prec.parameters["nonzero_initial_guess"] = False

###############################################################################
# 11. Solve the inverse problem using the Inexact (Gauss)-Newton CG method
#

# define parameters for the optimization
tol = 1e-7      # relative tolerance
c = 1e-5        # c_armijo for backtracking 
maxiter = 20    # Maximum number of (Gauss) Newton iterations
use_GaussNewton = False # If True we will use the Gauss-Newton Hessian,
                        # if False we will use the true Hessian

# initialize iter counters
iter = 1
total_cg_iter = 0
converged = False

# initialize the finite element vector to store the gradient and the Newton direction
g, a_delta = Vector(), Vector()
R.init_vector(a_delta,0)
R.init_vector(g,0)

du, dp = Vector(), Vector()
W.init_vector(du,1)
W.init_vector(dp,0)

a_prev = Function(Va)

print "Nit   CGit   cost          misfit        reg           sqrt(-G*D)    ||grad||       alpha  tolcg"

while iter <  maxiter and not converged:

    # assemble matrix C
    C =  assemble(C_equ)

    # solve the adoint problem
    adjoint_A, adjoint_RHS = assemble_system(a_adj, L_adj, bc_adj)
    solve(adjoint_A, p.vector(), adjoint_RHS)

    # assemble W_ua and R
    Wua = assemble (Wua_equ)
    Raa = assemble (Raa_equ)

    # evaluate the  gradient
    CT_p = Vector()
    C.init_vector(CT_p,1)
    C.transpmult(p.vector(), CT_p)
    MG = CT_p + R * a.vector()
    solve(M, g, MG)

    # calculate the norm of the gradient
    grad2 = g.inner(MG)
    gradnorm = sqrt(grad2)

    # set the CG tolerance (use Eisenstat - Walker termination criterion)
    if iter == 1:
        gradnorm_ini = gradnorm
    tolcg = min(0.5, sqrt(gradnorm/gradnorm_ini))

    # define the Hessian apply operator (with preconditioner)
    if use_GaussNewton:
        Hess_Apply = GaussNewtonHessian(R, C, state_A, adjoint_A, W )
        solver = CGSolver()
    else:
        Hess_Apply = Hessian(R, Raa, C, state_A, adjoint_A, W, Wua )
        solver = CGSolverSteihaug()
    
    solver.set_operator(Hess_Apply)
    solver.set_preconditioner(prec)
    # Set the tolerance for CG determined by the  Eisenstat-Walker conditions
    solver.parameters["rel_tolerance"] = tolcg

    # solve the Newton system H a_delta = - MG
    solver.solve(a_delta, -MG)
    total_cg_iter += solver.iter
    
    # linesearch
    alpha = 1
    descent = False
    no_backtrack = 0
    a_prev.assign(a)
    while not descent and no_backtrack < 10:
        a.vector().axpy(alpha, a_delta )

        # solve the state/forward problem
        state_A, state_b = assemble_system(a_state, L_state, bc_state)
        solve(state_A, u.vector(), state_b)

        # evaluate cost
        [cost_new, misfit_new, reg_new] = cost(u, ud, a, W, R)

        # check if Armijo conditions are satisfied
        if cost_new < cost_old + alpha * c * MG.inner(a_delta):
            cost_old = cost_new
            descent = True
        else:
            no_backtrack += 1
            alpha *= 0.5
            a.assign(a_prev)  # reset a
            
    if not descent:
        raise NameError("Backtracking failed")

    # calculate sqrt(-G * D)
    graddir = sqrt(- MG.inner(a_delta) )

    sp = ""
    print "%2d %2s %2d %3s %8.5e %1s %8.5e %1s %8.5e %1s %8.5e %1s %8.5e %1s %5.2f %1s %5.3e" % \
        (iter, sp, solver.iter, sp, cost_new, sp, misfit_new, sp, reg_new, sp, \
         graddir, sp, gradnorm, sp, alpha, sp, tolcg)
    
    # check for convergence
    if gradnorm < tol*gradnorm_ini and iter > 1:
        converged = True
        if use_GaussNewton:
            print "Gauss Newton's method converged in ",iter,"  iterations"
        else:
            print "Newton's method converged in ",iter,"  iterations"
        print "Total number of CG iterations: ", total_cg_iter
        
    iter += 1
    
if not converged:
    print "Newton's method did not converge in ", maxiter, " iterations"




print "Time elapsed: ", time.clock()-start

nb.multi1_plot([atrue, a], ["atrue", "a"])
nb.multi1_plot([u,p], ["u","p"], same_colorbar=False)

plt.show()

