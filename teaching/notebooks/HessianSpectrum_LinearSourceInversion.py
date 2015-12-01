from dolfin import *

import numpy as np
import time
import logging

import matplotlib.pyplot as plt
import nb

import sys
sys.path.append("../../")
from hippylib import *

start = time.clock()

logging.getLogger('FFC').setLevel(logging.WARNING)
logging.getLogger('UFL').setLevel(logging.WARNING)
set_log_active(False)

def pde_varf(u,a,p):
    return k*inner(nabla_grad(u), nabla_grad(p))*dx \
           + inner(nabla_grad(u), v*p)*dx \
           + c*u*p*dx \
           - a*p*dx

def u_boundary(x, on_boundary):
    return on_boundary and x[1] < DOLFIN_EPS

def solve(nx,ny, targets, rel_noise, gamma, delta, verbose=True):
    np.random.seed(seed=2)
    mesh = UnitSquareMesh(nx, ny)
    Vh1 = FunctionSpace(mesh, 'Lagrange', 1)
    
    Vh = [Vh1, Vh1, Vh1]
    if verbose:
        print "Number of dofs: STATE={0}, PARAMETER={1}, ADJOINT={2}".format(Vh[STATE].dim(), Vh[PARAMETER].dim(), Vh[ADJOINT].dim())


    u_bdr = Expression("0.0")
    u_bdr0 = Expression("0.0")
    bc = DirichletBC(Vh[STATE], u_bdr, u_boundary)
    bc0 = DirichletBC(Vh[STATE], u_bdr0, u_boundary)

    atrue = interpolate( Expression("exp( -50*(x[0] - .5)*(x[0] - .5) - 50*(x[1] - .5)*(x[1] - .5))"), Vh[PARAMETER]).vector()
    a0 = interpolate(Expression("0.0"), Vh[PARAMETER]).vector()
    
    pde = PDEVariationalProblem(Vh, pde_varf, bc, bc0)
 
    if verbose:
        print "Number of observation points: {0}".format(targets.shape[0])
        
    misfit = PointwiseStateObservation(Vh[STATE], targets)
    
    reg = LaplacianPrior(Vh[PARAMETER], gamma, delta)
                    
    #Generate synthetic observations
    utrue = pde.generate_state()
    x = [utrue, atrue, None]
    pde.solveFwd(x[STATE], x, 1e-9)
    misfit.B.mult(x[STATE], misfit.d)
    MAX = misfit.d.norm("linf")
    noise_std_dev = rel_noise * MAX
    randn_perturb(misfit.d, noise_std_dev)
    misfit.noise_variance = noise_std_dev*noise_std_dev

    if verbose:
        plt.figure(figsize=(18,4))
        nb.plot(Function(Vh[PARAMETER], atrue), mytitle = "True source", subplot_loc=131)
        nb.plot(Function(Vh[STATE], utrue), mytitle="True state", subplot_loc=132)
        nb.plot_pts(targets, misfit.d,mytitle="Observations", subplot_loc=133)
        plt.show()
    
    model = Model(pde, reg, misfit)
    u = model.generate_vector(STATE)
    a = a0.copy()
    p = model.generate_vector(ADJOINT)
    x = [u,a,p]
    mg = model.generate_vector(PARAMETER)
    model.solveFwd(u, x)
    model.solveAdj(p, x)
    model.evalGradientParameter(x, mg)
    model.setPointForHessianEvaluations(x)

    H = ReducedHessian(model, 1e-12)

    solver = CGSolverSteihaug()
    solver.set_operator(H)
    solver.set_preconditioner( reg.Rsolver )
    solver.parameters["print_level"] = -1
    solver.parameters["rel_tolerance"] = 1e-9
    solver.solve(a, -mg)

    if solver.converged:
        if verbose:
            print "CG converged in ", solver.iter, " iterations."
    else:
        print "CG did not converged."
        raise

    model.solveFwd(u, x, 1e-12)
 
    total_cost, reg_cost, misfit_cost = model.cost(x)

    if verbose:
        plt.figure(figsize=(18,4))
        nb.plot(Function(Vh[PARAMETER], a), mytitle = "Reconstructed source", subplot_loc=131)
        nb.plot(Function(Vh[STATE], u), mytitle="Reconstructed state", subplot_loc=132)
        nb.plot_pts(targets, misfit.B*u - misfit.d, mytitle="Misfit", subplot_loc=133)
        plt.show()

    H.misfit_only = True
    k_evec = 80
    p_evec = 5
    if verbose:
        print "Double Pass Algorithm. Requested eigenvectors: {0}; Oversampling {1}.".format(k_evec,p_evec)
    Omega = np.random.randn(a.array().shape[0], k_evec+p_evec)
    d, U = doublePassG(H, reg.R, reg.Rsolver, Omega, k_evec)

    if verbose:
        plt.figure()
        nb.plot_eigenvalues(d, mytitle="Generalized Eigenvalues")
        nb.plot_eigenvectors(Vh[PARAMETER], U, mytitle="Eigenvectors", which=[0,1,2,5,10,15])
        plt.show()
        
    return d, U, Vh[PARAMETER], solver.iter


ndim = 2
nx = 32
ny = 32

ntargets = 300
np.random.seed(seed=1)
targets = np.random.uniform(0.1,0.9, [ntargets, ndim] )
rel_noise = 0.01

gamma = 70.
delta = 1e-1

k = Expression("1.0")
v = Expression(("0.0", "0.0"))
c = Expression("0.")

d, U, Va, nit = solve(nx,ny, targets, rel_noise, gamma, delta)

gamma = 70.
delta = 1e-1

k = Expression("1.0")
v = Expression(("0.0", "0.0"))
c = Expression("0.")

n = [16,32,64]
d1, U1, Va1, niter1 = solve(n[0],n[0], targets, rel_noise, gamma, delta,verbose=False)
d2, U2, Va2, niter2 = solve(n[1],n[1], targets, rel_noise, gamma, delta,verbose=False)
d3, U3, Va3, niter3 = solve(n[2],n[2], targets, rel_noise, gamma, delta,verbose=False)

print "Number of Iterations: ", niter1, niter2, niter3
plt.figure(figsize=(18,4))
nb.plot_eigenvalues(d1, mytitle="Eigenvalues Mesh {0} by {1}".format(n[0],n[0]), subplot_loc=131)
nb.plot_eigenvalues(d2, mytitle="Eigenvalues Mesh {0} by {1}".format(n[1],n[1]), subplot_loc=132)
nb.plot_eigenvalues(d3, mytitle="Eigenvalues Mesh {0} by {1}".format(n[2],n[2]), subplot_loc=133)

nb.plot_eigenvectors(Va1, U1, mytitle="Mesh {0} by {1} Eigen".format(n[0],n[0]), which=[0,1,5])
nb.plot_eigenvectors(Va2, U2, mytitle="Mesh {0} by {1} Eigen".format(n[1],n[1]), which=[0,1,5])
nb.plot_eigenvectors(Va3, U3, mytitle="Mesh {0} by {1} Eigen".format(n[2],n[2]), which=[0,1,5])

plt.show()

gamma = 70.
delta = 1e-1

k = Expression("1.0")
v = Expression(("0.0", "0.0"))
c = Expression("0.")

rel_noise = [1e-3,1e-2,1e-1]
d1, U1, Va1, niter1 = solve(nx,ny, targets, rel_noise[0], gamma, delta,verbose=False)
d2, U2, Va2, niter2 = solve(nx,ny, targets, rel_noise[1], gamma, delta,verbose=False)
d3, U3, Va3, niter3 = solve(nx,ny, targets, rel_noise[2], gamma, delta,verbose=False)

print "Number of Iterations: ", niter1, niter2, niter3
plt.figure(figsize=(18,4))
nb.plot_eigenvalues(d1, mytitle="Eigenvalues rel_noise {0:g}".format(rel_noise[0]), subplot_loc=131)
nb.plot_eigenvalues(d2, mytitle="Eigenvalues rel_noise {0:g}".format(rel_noise[1]), subplot_loc=132)
nb.plot_eigenvalues(d3, mytitle="Eigenvalues rel_noise {0:g}".format(rel_noise[2]), subplot_loc=133)

nb.plot_eigenvectors(Va1, U1, mytitle="rel_noise {0:g} Eigen".format(rel_noise[0]), which=[0,1,5])
nb.plot_eigenvectors(Va2, U2, mytitle="rel_noise {0:g} Eigen".format(rel_noise[1]), which=[0,1,5])
nb.plot_eigenvectors(Va3, U3, mytitle="rel_noise {0:g} Eigen".format(rel_noise[2]), which=[0,1,5])

plt.show()

rel_noise = 0.01

k = Expression("1.0")
v = Expression(("0.0", "0.0"))
c = Expression("1.0")

d1, U1, Va1, niter1 = solve(nx,ny, targets, rel_noise, gamma, delta,verbose=False)
k = Expression("0.1")
d2, U2, Va2, niter2 = solve(nx,ny, targets, rel_noise, gamma, delta,verbose=False)
k = Expression("0.01")
d3, U3, Va3, niter3 = solve(nx,ny, targets, rel_noise, gamma, delta,verbose=False)

print "Number of Iterations: ", niter1, niter2, niter3
plt.figure(figsize=(18,4))
nb.plot_eigenvalues(d1, mytitle="Eigenvalues k=1.0", subplot_loc=131)
nb.plot_eigenvalues(d2, mytitle="Eigenvalues k=0.1", subplot_loc=132)
nb.plot_eigenvalues(d3, mytitle="Eigenvalues k=0.01", subplot_loc=133)

nb.plot_eigenvectors(Va1, U1, mytitle="k=1. Eigen", which=[0,1,5])
nb.plot_eigenvectors(Va2, U2, mytitle="k=0.1 Eigen", which=[0,1,5])
nb.plot_eigenvectors(Va3, U3, mytitle="k=0.01 Eigen", which=[0,1,5])

plt.show()
