import dolfin as dl
import sys
sys.path.append( "../../" )
from hippylib import *
import numpy as np
import matplotlib.pyplot as plt
sys.path.append("../poisson/")
from model_subsurf import Poisson, true_model
from posteriorDistribution import PosteriorDistribution, GaussianDistribution

if __name__ == "__main__":
    dl.set_log_active(False)
    sep = "\n"+"#"*80+"\n"
    print sep, "Set up the mesh and finite element spaces", sep
    ndim = 2
    nx = 64
    ny = 64
    mesh = dl.UnitSquareMesh(nx, ny)
    Vh2 = dl.FunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
    Vh = [Vh2, Vh1, Vh2]
    print "Number of dofs: STATE={0}, PARAMETER={1}, ADJOINT={2}".format(Vh[STATE].dim(), Vh[PARAMETER].dim(), Vh[ADJOINT].dim())
    
    print sep, "Set up the location of observation, Prior Information, and model", sep
    ntargets = 5
    np.random.seed(seed=1)
    targets = np.random.uniform(0.1,0.9, [ntargets, ndim] )
    print "Number of observation points: {0}".format(ntargets)
    
    gamma = .1
    delta = .5
    
    anis_diff = dl.Expression(code_AnisTensor2D)
    anis_diff.theta0 = 2.
    anis_diff.theta1 = .5
    anis_diff.alpha = math.pi/4
    atrue = true_model(Vh[PARAMETER], gamma, delta,anis_diff)
        
    locations = np.array([[0.1, 0.1], [0.1, 0.9], [.5,.5], [.9, .1], [.9, .9]])
    if 1:
        pen = 0.001
        prior = MollifiedBiLaplacianPrior(Vh[PARAMETER], gamma, delta, locations, atrue, anis_diff, pen)
    else:
        pen = 1e4
        prior = ConstrainedBiLaplacianPrior(Vh[PARAMETER], gamma, delta, locations, atrue, anis_diff, pen)
        
    print "Prior regularization: (delta_x - gamma*Laplacian)^order: delta={0}, gamma={1}, order={2}".format(delta, gamma,2)    
            
    model = Poisson(mesh, Vh, targets, prior)
    
    #Generate synthetic observations
    utrue = model.generate_vector(STATE)
    x = [utrue, atrue, None]
    model.solveFwd(x[STATE], x, 1e-9)
    model.B.mult(x[STATE], model.u_o)
    rel_noise = 0.2
    MAX = model.u_o.norm("linf")
    noise_std_dev = rel_noise * MAX
    randn_perturb(model.u_o, noise_std_dev)
    model.noise_variance = noise_std_dev*noise_std_dev
   
    print sep, "Test the gradient and the Hessian of the model", sep
    a0 = dl.interpolate(dl.Expression("sin(x[0])"), Vh[PARAMETER])
    modelVerify(model, a0.vector(), 1e-4, 1e-6)

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
    plt.figure()
    plt.semilogy(d)
    posterior.mean = x[PARAMETER]
    
        
    noise = model.generate_vector(PARAMETER)
    noise_size = noise.array().shape[0]
    
    start = posterior.mean.copy()
    
    dir = model.generate_vector(PARAMETER)
    noise.set_local( np.random.randn( noise_size ) )
    posterior.sample(noise, dir, add_mean=False)
    
    pdf = PosteriorDistribution(model, solver.final_cost)
    pdfg = GaussianDistribution(posterior.Hlr, posterior.mean)
    pdfp = GaussianDistribution(prior.R, posterior.mean)
    
    check_derivatives = False
    if check_derivatives:
        n = 20
        all_eps = 1e-3*np.power(2., -np.arange(0,n))
        err_grad = np.ndarray(all_eps.shape)
        err_H = np.ndarray(all_eps.shape)
        grad = model.generate_vector(PARAMETER)
        pdf.gradient(start, grad)
        graddir = grad.inner(dir)
        c = pdf(start)
        for i in range(n):
            eps = all_eps[i]
            aplus = start.copy()
            aplus.axpy(eps, dir)        
            cplus = pdf(aplus)
        
            err_grad[i] =  np.abs( (cplus-c)/eps - graddir)

    
        Hdir = model.generate_vector(PARAMETER)
        pdf.hessian_apply(start, dir, Hdir)
        grad_plus = model.generate_vector(PARAMETER)
        for i in range(n):
            eps = all_eps[i]
            aplus = start.copy()
            aplus.axpy(eps, dir)        
            pdf.gradient(aplus, grad_plus)        
            err = grad_plus.copy()
            err.axpy(-1, grad)
            err *= 1/eps
            err.axpy(-1, Hdir)
            err_H[i] = err.norm("linf")
        
        err_gradg = np.ndarray(all_eps.shape)
        err_Hg = np.ndarray(all_eps.shape)
        pdfg.gradient(start, grad)
        graddir = grad.inner(dir)
        c = pdf(start)
        for i in range(n):
            eps = all_eps[i]
            aplus = start.copy()
            aplus.axpy(eps, dir)        
            cplus = pdfg(aplus)
        
            err_gradg[i] =  np.abs( (cplus-c)/eps - graddir)

        pdfg.hessian_apply(start, dir, Hdir)
        grad_plus = model.generate_vector(PARAMETER)
        for i in range(n):
            eps = all_eps[i]
            aplus = start.copy()
            aplus.axpy(eps, dir)        
            pdfg.gradient(aplus, grad_plus)        
            err = grad_plus.copy()
            err.axpy(-1, grad)
            err *= 1/eps
            err.axpy(-1, Hdir)
            err_Hg[i] = err.norm("linf")
        
        plt.figure()
        plt.subplot(221)
        plt.loglog(all_eps, err_grad)
        plt.ylabel("Err Grad Post")
        plt.subplot(222)
        plt.loglog(all_eps, err_H)
        plt.ylabel("Err H Post")
        plt.subplot(223)
        plt.loglog(all_eps, err_gradg)
        plt.ylabel("Err Grad Gaussian")
        plt.subplot(224)
        plt.loglog(all_eps, err_Hg)
        plt.ylabel("Err H Gaussian")
    
    iact = np.where(d > .1)
    alpha = np.zeros(d.shape)
    alpha[iact] = np.random.randn( iact[0].shape[0] )
    print alpha
    dir = model.generate_vector(PARAMETER)
    dir.set_local( np.dot(U, alpha) )
#    l = 1/math.sqrt(1+d[i])
#    dir = model.generate_vector(PARAMETER)
#    noise.set_local( np.random.randn( noise_size ) )
#    posterior.sample(noise, dir, add_mean=False)
    H = ReducedHessian(model, solver.parameters["inner_rel_tolerance"], gauss_newton_approx=False, misfit_only=False)
    l = 1/math.sqrt(H.inner(dir, dir));
    
    npoints = 100
    d = np.linspace(-3*l, 3*l, npoints, endpoint=True)
    ct = np.ndarray(npoints)
    cg = np.ndarray(npoints)
    for i in range(npoints):
        pp = posterior.mean.copy()
        pp.axpy(d[i], dir)
        ct[i] = pdf(pp)
        cg[i] = pdfg(pp)

        
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(d, ct, "-b")
    plt.plot(d, cg, "--r")
    plt.subplot(1,2,2)
    plt.plot(d, ct-cg, "-b")
    plt.show()

