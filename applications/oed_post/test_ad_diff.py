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

import sys
sys.path.append( "../../" )
from hippylib import *
sys.path.append("../ad_diff/")
from model_gls import TimeDependentAD, computeVelocityField
from posteriorDistribution import PosteriorDistribution, GaussianDistribution

if __name__ == "__main__":
    dl.set_log_active(False)
    np.random.seed(1)
    sep = "\n"+"#"*80+"\n"
    print sep, "Set up the mesh and finite element spaces.\n","Compute wind velocity", sep
    mesh = dl.refine( dl.Mesh("../ad_diff/ad_20.xml") )
    wind_velocity = computeVelocityField(mesh)
    Vh = dl.FunctionSpace(mesh, "Lagrange", 2)
    print "Number of dofs: {0}".format( Vh.dim() )
    
    print sep, "Set up Prior Information and model", sep
    
    true_initial_condition = dl.interpolate(dl.Expression('min(0.5,exp(-100*(pow(x[0]-0.35,2) +  pow(x[1]-0.7,2))))'), Vh).vector()

    orderPrior = 2
    
    if orderPrior == 1:
        gamma = 1
        delta = 1e1
        prior = LaplacianPrior(Vh, gamma, delta)
    elif orderPrior == 2:
        gamma = 1
        delta = 8
        prior = BiLaplacianPrior(Vh, gamma, delta)
        
#    prior.mean = interpolate(Expression('min(0.6,exp(-50*(pow(x[0]-0.34,2) +  pow(x[1]-0.71,2))))'), Vh).vector()
    prior.mean = dl.interpolate(dl.Expression('0.5'), Vh).vector()
    
    print "Prior regularization: (delta - gamma*Laplacian)^order: delta={0}, gamma={1}, order={2}".format(delta, gamma,orderPrior)

    problem = TimeDependentAD(mesh, [Vh,Vh,Vh], 0., 4., 1., .2, wind_velocity, True, prior)
    
    print sep, "Generate synthetic observation", sep
    rel_noise = 0.001
    utrue = problem.generate_vector(STATE)
    x = [utrue, true_initial_condition, None]
    problem.solveFwd(x[STATE], x, 1e-9)
    MAX = utrue.norm("linf", "linf")
    noise_std_dev = rel_noise * MAX
    problem.ud.copy(utrue)
    problem.ud.randn_perturb(noise_std_dev)
    problem.noise_variance = noise_std_dev*noise_std_dev
    
    print sep, "Test the gradient and the Hessian of the model", sep
    a0 = true_initial_condition.copy()
    modelVerify(problem, a0, 1e-4, 1e-4)
    
    print sep, "Compute the reduced gradient and hessian", sep
    [u,a,p] = problem.generate_vector()
    problem.solveFwd(u, [u,a,p], 1e-12)
    problem.solveAdj(p, [u,a,p], 1e-12)
    mg = problem.generate_vector(PARAMETER)
    grad_norm = problem.evalGradientParameter([u,a,p], mg)
        
    print "(g,g) = ", grad_norm
    
    H = ReducedHessian(problem, 1e-12, gauss_newton_approx=False, misfit_only=True) 
    
    print sep, "Compute the low rank Gaussian Approximation of the posterior", sep   
    k = 80
    p = 20
    print "Double Pass Algorithm. Requested eigenvectors: {0}; Oversampling {1}.".format(k,p)
    Omega = np.random.randn(a.array().shape[0], k+p)
    d, U = singlePassG(H, prior.R, prior.Rsolver, Omega, k, check_Bortho=False, check_Aortho=False, check_residual=False)
    posterior = GaussianLRPosterior( prior, d, U )
    
    print sep, "Find the MAP point", sep
    
    H.misfit_only = False
        
    solver = CGSolverSteihaug()
    solver.set_operator(H)
    solver.set_preconditioner( posterior.Hlr )
    solver.parameters["print_level"] = 1
    solver.parameters["rel_tolerance"] = 1e-6
    solver.solve(a, -mg)
    problem.solveFwd(u, [u,a,p], 1e-12)
 
    total_cost, reg_cost, misfit_cost = problem.cost([u,a,p])
    print "Total cost {0:5g}; Reg Cost {1:5g}; Misfit {2:5g}".format(total_cost, reg_cost, misfit_cost)
    
    posterior.mean = a
    
    pdf = PosteriorDistribution(problem, total_cost)
    pdfg = GaussianDistribution(posterior.Hlr, posterior.mean)
    pdfp = GaussianDistribution(prior.R, posterior.mean)
    
    i = 0
    dir = problem.generate_vector(PARAMETER)
    dir.set_local(U[:,i])
    
    npoints = 1000
    d = np.linspace(0, 1, npoints, endpoint=True)
    ct = np.ndarray(npoints)
    cg = np.ndarray(npoints)
    cp = np.ndarray(npoints)
    for i in range(npoints):
        pp = posterior.mean.copy()
        pp.axpy(d[i], dir)
        ct[i] = pdf(pp)
        cg[i] = pdfg(pp)
        cp[i] = pdfp(pp)

        
    plt.figure()
    plt.plot(d, ct, "-b")
    plt.plot(d, cg, "-.g")
    plt.plot(d, cp, "--r")
    plt.show()
