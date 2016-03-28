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
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

import sys
sys.path.append( "../../" )
from hippylib import *

class MyOperator:
    def __init__(self, op,op2):
        self.op = op
        self.op2 = op2
        self.op2x = dl.Vector()
        self.z = dl.Vector()
        self.op.init_vector(self.op2x,1)
        self.op.init_vector(self.z,0)
            
    def init_vector(self,x,dim):
        self.op2.init_vector(x,1)
            
    def mult(self,x,y):
        self.op2.mult(x,self.op2x)
        self.op.solve(self.z,self.op2x)
        self.op2.mult(self.z,y)

def true_model(Vh, gamma, delta, anis_diff):
    prior = BiLaplacianPrior(Vh, gamma, delta, anis_diff )
    noise = dl.Vector()
    prior.init_vector(noise,"noise")
    noise_size = noise.array().shape[0]
    noise.set_local( np.random.randn( noise_size ) )
    atrue = dl.Vector()
    prior.init_vector(atrue, 0)
    prior.sample(noise,atrue)
    return atrue

def u_boundary(x, on_boundary):
    return on_boundary and ( x[1] < dl.DOLFIN_EPS or x[1] > 1.0 - dl.DOLFIN_EPS)

if __name__ == "__main__":
    dl.set_log_active(False)
    sep = "\n"+"#"*80+"\n"
    print sep, "Set up the mesh and finite element spaces", sep
    ndim = 2
    nx = 16
    ny = 16
    mesh = dl.UnitSquareMesh(nx, ny)
    Vh2 = dl.FunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
    Vh = [Vh2, Vh1, Vh2]
    print "Number of dofs: STATE={0}, PARAMETER={1}, ADJOINT={2}".format(Vh[STATE].dim(), Vh[PARAMETER].dim(), Vh[ADJOINT].dim())
    
    # Initialize Expressions
    f = dl.Expression("0.0")
        
    u_bdr = dl.Expression("x[1]")
    u_bdr0 = dl.Expression("0.0")
    bc = dl.DirichletBC(Vh[STATE], u_bdr, u_boundary)
    bc0 = dl.DirichletBC(Vh[STATE], u_bdr0, u_boundary)
    
    def pde_varf(u,a,p):
        return dl.exp(a)*(dl.Constant(1.)+u*u)*dl.inner(dl.nabla_grad(u), dl.nabla_grad(p))*dl.dx - f*p*dl.dx
    
    pde = PDEVariationalProblem(Vh, pde_varf, bc, bc0)
 
    ntargets = 10
    np.random.seed(seed=1)
    targets = np.random.uniform(0.1,0.9, [ntargets, ndim] )
    print "Number of observation points: {0}".format(ntargets)
    misfit = PointwiseStateObservation(Vh[STATE], targets)
    
    
    gamma = .1
    delta = .5
    
    anis_diff = dl.Expression(code_AnisTensor2D)
    anis_diff.theta0 = 2.
    anis_diff.theta1 = .5
    anis_diff.alpha = math.pi/4
    atrue = true_model(Vh[PARAMETER], gamma, delta,anis_diff)
        
    locations = np.array([[0.1, 0.1], [0.1, 0.9], [.5,.5], [.9, .1], [.9, .9]])

    pen = 1e-2
    prior = MollifiedBiLaplacianPrior(Vh[PARAMETER], gamma, delta, locations, atrue, anis_diff, pen)
        
    print "Prior regularization: (delta_x - gamma*Laplacian)^order: delta={0}, gamma={1}, order={2}".format(delta, gamma,2)    
                
    #Generate synthetic observations
    utrue = pde.generate_state()
    x = [utrue, atrue, None]
    pde.solveFwd(x[STATE], x, 1e-9)
    misfit.B.mult(x[STATE], misfit.d)
    rel_noise = 0.05
    MAX = misfit.d.norm("linf")
    noise_std_dev = rel_noise * MAX
    randn_perturb(misfit.d, noise_std_dev)
    misfit.noise_variance = noise_std_dev*noise_std_dev
    
    model = Model(pde,prior, misfit)
    
    print sep, "Find the MAP point", sep
    a0 = prior.mean.copy()
    solver = ReducedSpaceNewtonCG(model)
    solver.parameters["rel_tolerance"] = 1e-9
    solver.parameters["abs_tolerance"] = 1e-12
    solver.parameters["max_iter"]      = 25
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["c_armijo"] = 1e-4
    solver.parameters["GN_iter"] = 5
    
    x = solver.solve(a0)
    
    if solver.converged:
        print "\nConverged in ", solver.it, " iterations."
    else:
        print "\nNot Converged"

    print "Termination reason: ", solver.termination_reasons[solver.reason]
    print "Final gradient norm: ", solver.final_grad_norm
    print "Final cost: ", solver.final_cost
        
    print sep, "Compute the low rank Gaussian Approximation of the posterior", sep
    model.setPointForHessianEvaluations(x)
    Hmisfit = ReducedHessian(model, solver.parameters["inner_rel_tolerance"], gauss_newton_approx=False, misfit_only=True)
    k = min(50, ntargets)
    p = 20
    print "Double Pass Algorithm. Requested eigenvectors: {0}; Oversampling {1}.".format(k,p)
    Omega = np.random.randn(x[PARAMETER].array().shape[0], k+p)
    d, U = doublePassG(Hmisfit, prior.R, prior.Rsolver, Omega, k, check_Bortho=False, check_Aortho=False, check_residual=False)
    d[d < 0 ] = 0.
    posterior = GaussianLRPosterior(prior, d, U)
    posterior.mean = x[PARAMETER].copy()
    
    
    nsamples = 500
    my_u = model.generate_vector(STATE)
        
    print sep, "Generate samples from Prior and Posterior\n", sep
    fid_prior = dl.File("samples/sample_prior.pvd")
    fid_post  = dl.File("samples/sample_post.pvd")
    cost = np.zeros([nsamples,3])
    noise = dl.Vector()
    posterior.init_vector(noise,"noise")
    noise_size = noise.array().shape[0]
    s_prior = dl.Function(Vh[PARAMETER], name="sample_prior")
    s_post = dl.Function(Vh[PARAMETER], name="sample_post")
    for i in range(nsamples):
        noise.set_local( np.random.randn( noise_size ) )
        posterior.sample(noise, s_prior.vector(), s_post.vector())
        fid_prior << s_prior
        fid_post << s_post
        model.solveFwd(my_u, [my_u, s_post.vector()])
        cost[i,0], cost[i,1], cost[i,2] = model.cost([my_u, s_post.vector()])
    
    Gamma_post = to_dense( MyOperator(posterior.Hlr, prior.M) )
    Gamma_prior = to_dense( MyOperator(prior.Rsolver, prior.M) )
    M_dense = to_dense(prior.M)
    d_gaussianPost, U_gaussianPost = eigh(Gamma_post,M_dense)
    d_gaussianPost = d_gaussianPost[::-1]
    U_gaussianPost = U_gaussianPost[:,::-1]
    
    d_prior, U_prior = eigh(Gamma_prior,M_dense)
    d_prior = d_prior[::-1]
    U_prior = U_prior[:,::-1]
    
    print sep, "Generate samples from Prior and Posterior using KLE\n", sep
    fid_prior = dl.File("samplesKLE/sample_prior.pvd")
    fid_post  = dl.File("samplesKLE/sample_post.pvd")
    noise_size = d_gaussianPost.shape[0]
    s_prior = dl.Function(Vh[PARAMETER], name="sample_prior")
    s_post = dl.Function(Vh[PARAMETER], name="sample_post")
    costKLE = np.zeros([nsamples,3])
    for i in range(nsamples):
        eta = np.random.randn(noise_size)
        s_prior_data = U_prior.dot(eta*np.sqrt(d_prior))
        s_post_data  = U_gaussianPost.dot(eta*np.sqrt(d_gaussianPost))
        s_prior.vector().set_local(s_prior_data)
        s_post.vector().set_local(s_post_data)
        fid_prior << s_prior
        fid_post << s_post
        model.solveFwd(my_u, [my_u, s_post.vector()])
        costKLE[i,0], costKLE[i,1], costKLE[i,2] = model.cost([my_u, s_post.vector()])
        
    plt.figure()
    plt.subplot(121)
    plt.plot(cost)
    plt.subplot(122)
    plt.plot(costKLE)
    plt.show()


