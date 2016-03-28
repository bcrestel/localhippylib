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

import dolfin as dl
import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

import sys
sys.path.append( "../../" )
from hippylib import *
from posteriorDistribution import PosteriorDistribution, GaussianDistribution, check_derivatives_pdf

def exportU(U, Vh, fname, varname = "evect", normalize=1):

        evect = dl.Function(Vh, name=varname)
        fid = dl.File(fname)
        
        for i in range(0,U.shape[1]):
            Ui = U[:,i]
            if normalize:
                s = 1/np.linalg.norm(Ui, np.inf)
                evect.vector().set_local(s*Ui)
            else:
                evect.vector().set_local(Ui)
            fid << evect

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

def marginal_distribution(map, pdf, pdfg, l, dir):   
    npoints = 50
    r = np.linspace(-3*l, 3*l, npoints, endpoint=True)
    ct = np.ndarray(npoints)
    cg = np.ndarray(npoints)
    for i in range(npoints):
        pp = map.copy()
        pp.axpy(r[i], dir)
        ct[i] = pdf(pp)
        cg[i] = pdfg(pp)
            
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(r, ct, "-b", r, cg, "--r")
    plt.title( "Marginal distribution lambda = {0:f}".format(math.pow(l, -2.)-1) )
    plt.subplot(1,2,2)
    plt.plot(r, ct-cg, "-b")
    plt.show()

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
           
    print sep, "Test the gradient and the Hessian of the model", sep
    a0 = dl.interpolate(dl.Expression("sin(x[0])"), Vh[PARAMETER])
    if 0:
        modelVerify(model, a0.vector(), 1e-12)
        plt.show()

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
    #d, U = singlePassG(Hmisfit, model.R, model.Rsolver, Omega, k, check_Bortho=True, check_Aortho=True, check_residual=True)
    d, U = doublePassG(Hmisfit, prior.R, prior.Rsolver, Omega, k, check_Bortho=False, check_Aortho=False, check_residual=False)
    posterior = GaussianLRPosterior(prior, d, U)
    
    if False:
        plt.figure()
        plt.semilogy(d, "ob")
        plt.semilogy(np.ones(d.shape), "-r")
        plt.title("Generalized Eigenvalues")
        plt.show()
                
    posterior.mean = x[PARAMETER]
        
    start = posterior.mean.copy()
        
    pdf = PosteriorDistribution(model, solver.final_cost)
    pdfg = GaussianDistribution(posterior.Hlr, posterior.mean)
    pdfp = GaussianDistribution(prior.R, posterior.mean)
    
    check_derivatives = False
    if check_derivatives:
        noise = pdf.generate_vector()
        noise_size = noise.array().shape[0]
        dir = pdf.generate_vector()
        noise.set_local( np.random.randn( noise_size ) )
        posterior.sample(noise, dir, add_mean=False)
        check_derivatives_pdf(pdf,start, dir)
        check_derivatives_pdf(pdfg,start, dir)
    
    check_along_eigenfunctions = True
    check_random_eigen_combination = False
    check_random_dir = False
    
    if check_along_eigenfunctions:
        for i in range(k):
            dir = model.generate_vector(PARAMETER)
            dir.set_local(U[:,i])
            l = 1./math.sqrt(1.+d[i])
            marginal_distribution(posterior.mean, pdf, pdfg, l, dir)
            
    if check_random_eigen_combination:
        iact = np.where(d > .1)
        alpha = np.zeros(d.shape)
        alpha[iact] = np.random.randn( iact[0].shape[0] )
        alpha *= 1./np.linalg.norm(alpha)
        dir = model.generate_vector(PARAMETER)
        dir.set_local( np.dot(U, alpha) )
        H = ReducedHessian(model, solver.parameters["inner_rel_tolerance"], gauss_newton_approx=False, misfit_only=False)
        l = 1./math.sqrt(H.inner(dir, dir))
        marginal_distribution(posterior.mean, pdf, pdfg, l, dir)
        
    if check_random_dir:
        noise = model.generate_vector(PARAMETER)
        noise_size = noise.array().shape[0]
        noise.set_local( np.random.randn( noise_size ) )
        dir = model.generate_vector(PARAMETER)
        posterior.sample(noise, dir, add_mean=False)
        H = ReducedHessian(model, solver.parameters["inner_rel_tolerance"], gauss_newton_approx=False, misfit_only=False)
        l = 1./math.sqrt(H.inner(dir, dir))
        marginal_distribution(posterior.mean, pdf, pdfg, l, dir)
    
    class MyOperator:
        def __init__(self, op,op2):
            self.op = op
            self.op2 = op2
            
        def init_vector(self,x,dim):
            self.op2.init_vector(x,1)
            
        def mult(self,x,y):
            op2x = dl.Vector()
            z = dl.Vector()
            self.op2.init_vector(op2x,0)
            self.op.init_vector(z,0)
            self.op2.mult(x,op2x)
            self.op.solve(z,op2x)
            self.op2.mult(z,y)
    
    Gamma_post = to_dense( MyOperator(posterior.Hlr, prior.M) )
    M_dense = to_dense( prior.M )
    d1, U1 = eigh(Gamma_post,M_dense)
    d1 = d1[::-1]
    U1 = U1[:,::-1]
    exportU(U1, Vh[PARAMETER], "EPost/Eigen_Post.pvd")
        
    plt.figure()
    plt.semilogy(d1, "ob")
    plt.semilogy(np.ones(d1.shape), "-r")
    plt.title("Posterior Eigenvalues")
    plt.show()
    
    
    for i in range(k):
        dir = model.generate_vector(PARAMETER)
        dir.set_local(U1[:,i])
        l = 1./math.sqrt(posterior.Hlr.inner(dir,dir))
        marginal_distribution(posterior.mean, pdf, pdfg, l, dir)
    