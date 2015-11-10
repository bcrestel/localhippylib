## Steepest descent solution of a coefficient field inversion in an
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

# Import dependencies
from dolfin import *

import numpy as np
import time
import logging

start = time.clock()

logging.getLogger('FFC').setLevel(logging.WARNING)
logging.getLogger('UFL').setLevel(logging.WARNING)
set_log_active(False)

np.random.seed(seed=1)

# The cost function evaluation:
def cost(u, ud, a, W, R):
    diff = u.vector() - ud.vector()
    reg = 0.5 * a.vector().inner(R*a.vector() )
    misfit = 0.5 * diff.inner(W * diff)
    return [reg + misfit, misfit, reg]


# Model set up:

# create mesh and define function spaces
nx = 32
ny = 32
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, 'Lagrange', 1)
V2 = FunctionSpace(mesh, 'Lagrange', 2)

# The true and inverted parameter
atrue = interpolate(Expression('8. - 4.*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2) < pow(0.2,2))'), V)
a = interpolate(Expression("8."),V)

File("parameter_true.pvd") << atrue
File("parameter_initial_guess.pvd") << a

# define function for state and adjoint
u = Function(V2)
p = Function(V2)

# define Trial and Test Functions
u_trial, p_trial, a_trial = TrialFunction(V2), TrialFunction(V2), TrialFunction(V)
u_test, p_test, a_test = TestFunction(V2), TestFunction(V2), TestFunction(V)

# initialize input functions
f = Constant("1.0")
u0 = Constant("0.0")

# noise level
noise_level = 0.05

# define parameters for the optimization
tol = 1e-7
gamma = 1e-16
maxiter = 1000
plot_any = 30

# initialize iter counters
iter = 1
converged = False

# set up dirichlet boundary conditions
def u0_boundary(x,on_boundary):
    return on_boundary
bc2 = DirichletBC(V2, u0, u0_boundary)

# Set up synthetic observations:
# weak form for setting up the synthetic observations
a_goal = inner( atrue * nabla_grad(u_trial), nabla_grad(u_test)) * dx
L_goal = f * u_test * dx

# solve the forward/state problem to generate synthetic observations
goal_A, goal_b = assemble_system(a_goal, L_goal, bc2)

utrue = Function(V2)
solve(goal_A, utrue.vector(), goal_b)

ud = Function(V2)
ud.assign(utrue)

# perturb state solution and create synthetic measurements ud
# ud = u + ||u||/SNR * random.normal
MAX = ud.vector().norm("linf")
noise = Vector()
goal_A.init_vector(noise,1)
noise.set_local( noise_level * MAX * np.random.normal(0, 1, len(ud.vector().array())) )
bc2.apply(noise)

ud.vector().axpy(1., noise)

File("state_true.pvd") << utrue
File("observation.pvd") << ud

# Setting up the state equations, right hand side for the adjoint and
# the neccessary matrices:

# weak form for setting up the state equation
a_state = inner( a * nabla_grad(u_trial), nabla_grad(u_test)) * dx
L_state = f * u_test * dx
W_equ   = inner(u_trial, u_test) * dx

# weak form for setting up the right hand side of the adjoint
L_adjoint = -inner(u - ud, u_test) * dx

# weak form for setting up matrices
CT_equ  = inner(a_test * nabla_grad(u), nabla_grad(p_trial)) * dx
M_equ   = inner(a_trial, a_test) * dx
R_equ   = gamma * inner(nabla_grad(a_trial), nabla_grad(a_test)) * dx

# assemble matrices M, W, and R
M = assemble(M_equ)
W = assemble(W_equ)
R = assemble(R_equ)

# solve state equation
A, state_b = assemble_system (a_state, L_state, bc2)
solve (A, u.vector(), state_b)

# evaluate cost
[cost_old, misfit_old, reg_old] = cost(u, ud, a, W, R)

# The steepest descent with Armijo line search:

# initializations
g = Vector()
R.init_vector(g,0)

a_prev = Function(V)

print "Nit  cost          misfit        reg         ||grad||       alpha  "

while iter <  maxiter and not converged:

    # assemble matrix C
    CT =  assemble(CT_equ)

    # solve the adoint problem
    A, adjoint_RHS = assemble_system(a_state, L_adjoint, bc2)
    solve(A, p.vector(), adjoint_RHS)

    # evaluate the  gradient
    MG = CT*p.vector() + R * a.vector()
    solve(M, g, MG)

    # calculate the norm of the gradient
    grad_norm2 = g.inner(MG)
    gradnorm = sqrt(grad_norm2)

    # linesearch
    c_armijo = 1e-4
    alpha = 1e5
    it_backtrack = 0
    a_prev.assign(a)
    for it_backtrack in range(20):

        a.vector().axpy(-alpha, g )

        # solve the state/forward problem
        state_A, state_b = assemble_system(a_state, L_state, bc2)
        solve(state_A, u.vector(), state_b)

        # evaluate cost
        [cost_new, misfit_new, reg_new] = cost(u, ud, a, W, R)

        # check if Armijo conditions are satisfied
        if cost_new < cost_old - alpha * c_armijo * grad_norm2:
            cost_old = cost_new
            break
        else:
            alpha *= 0.5
            a.assign(a_prev)  # reset a

    sp = ""
    print "%3d %1s %8.5e %1s %8.5e %1s %8.5e %1s %8.5e %1s %5.2f" % \
        (iter, sp, cost_new, sp, misfit_new, sp, reg_new, sp, \
        gradnorm, sp, alpha)

    # check for convergence
    if gradnorm < tol and iter > 1:
        converged = True
        print "Steepest descent converged in ",iter,"  iterations"

    iter += 1

if not converged:
    print "Steepest descent method did not converge in ", maxiter, " iterations"

print "Time elapsed: ", time.clock()-start

File("parameter_inverted.pvd") << a
File("state_inverted.pvd") << u
File("adjoint_inverted.pvd") << p
