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
import sys
sys.path.append( "../../" )
from hippylib import *
import numpy as np
import matplotlib.pyplot as plt

def true_model(Vh):
    cpp =
    """
    class BetaTrue : public Expression
    {
    public:

    BetaTrue() :
    Expression()
    {
    }

void eval(Array<double>& values, const Array<double>& x) const
  {
  double pi = 4.*atan(1.);
  double val = 0.;
  double x0 = x[0];
  if(x0 < 3750.)
     val = 1000.*(1. + sin(2*pi*x0/5000.));
  else if (x0 < 4000.)
     val = 1000.*(16 - x0/250.);
  else
     val = 1000.;
  values[0] = log(val);
  }  
};
    """
    return dl.interpolate(dl.Expression(cpp), Vh).vector()

class Arolla:
    def __init__(self, mesh, Vh, targets, prior):
        """
        Construct a model by proving
        - the mesh
        - the finite element spaces for the STATE/ADJOINT variable and the PARAMETER variable
        - the Prior information
        """
        self.mesh = mesh
        self.Vh = Vh
        
        self.rho = dl.Constant(910.) #kg/m^3
        self.g   = dl.Constant(9.81) #m/s^2
        self.n   = dl.Constant(3)
        self.RateFactor = dl.Constant(1e-16) #Pa^-n years^-1 
        
        # Initialize Expressions
        self.f_vel = dl.Expression([dl.Constant(0.), -self.rho*self.g])
                        
        # Assemble constant matrices      
        self.prior = prior
        self.B = assemblePointwiseObservation(Vh[STATE],targets)
                
        self.A = []
        self.At = []
        self.C = []
        self.Raa = []
        self.Wau = []
        
        self.u_o = dl.Vector()
        self.B.init_vector(self.u_o,0)
        self.noise_variance = 0
        
    def _forwardvarf(self, a, up, vq):
        n = dl.FacetNormal(self.mesh)
        h = dl.CellSize(self.mesh)
        pen = dl.Constant(10)
        
        uh,ph = dl.split(up)
        vh,qh = dl.split(vq)
        
        eta = dl.power(self.RateFactor, -1./self.n)*dl.power( dl.sym( dl.grad(uh) ) , .5*( dl.Constant(1.) - self.n)/self.n )
        def t(u,p,n): return dl.dot(eta*dl.sym(dl.grad(u)),n) - p*n
        
        a11 = eta*dl.inner(dl.sym( dl.grad(uh) ), dl.sym( dl.grad(vh) ))*dl.dx
        a12 = -dl.nabla_div(vh)*ph*dl.dx
        a21 = -dl.nabla_div(uh)*qh*dl.dx
        robinvarf  = dl.exp(a)*dl.inner(uh - dl.dot(dl.outer(n,n),uh), vh) * self.ds(1)
        weak_bc = -dl.dot(n, t(uh, ph, n) )*dl.dot(vh,n)*dl.ds(1) - dl.dot(n, t(vh, qh, n) )*dl.dot(uh,n)*dl.ds(1) + pen/h*dl.dot(uh,n)*dl.dot(vh, n)*dl.ds(1)
        
        return a11+a12+a21+robinvarf+weak_bc
        
    def generate_vector(self, component="ALL"):
        """
        Return the list x=[u,a,p] where:
        - u is any object that describes the state variable
        - a is a Vector object that describes the parameter variable.
          (Need to support linear algebra operations)
        - p is any object that describes the adjoint variable
        
        If component is STATE, PARAMETER, or ADJOINT return x[component]
        """
        if component == "ALL":
            x = [dl.Vector(), dl.Vector(), dl.Vector()]
            self.B.init_vector(x[STATE],1)
            self.prior.init_vector(x[PARAMETER],0)
            self.B.init_vector(x[ADJOINT], 1)
        elif component == STATE:
            x = dl.Vector()
            self.B.init_vector(x,1)
        elif component == PARAMETER:
            x = dl.Vector()
            self.prior.init_vector(x,0)
        elif component == ADJOINT:
            x = dl.Vector()
            self.B.init_vector(x,1)
            
        return x
    
    def init_parameter(self, a):
        """
        Reshape a so that it is compatible with the parameter variable
        """
        self.prior.init_vector(a,0)
        
    def assembleA(self,x, assemble_adjoint = False):
        """
        Assemble the matrices and rhs for the forward/adjoint problems
        """
        trial = dl.TrialFunction(self.Vh[STATE])
        test = dl.TestFunction(self.Vh[STATE])
        a = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        up = dl.Function(self.Vh[PARAMETER], x[STATE])
        res_varf = self._forwardvarf(self, a, up, test)
        Avarf = dl.derivative(res_varf, up, trial)
        if not assemble_adjoint:
            A = dl.assemble(Avarf)
        else:
            # Assemble the adjoint of A (i.e. the transpose of A)
            A = dl.assemble(dl.adjoint(Avarf))
            
        return A
    
    def assembleC(self, x):
        """
        Assemble the derivative of the forward problem with respect to the parameter
        """
        trial = dl.TrialFunction(self.Vh[PARAMETER])
        test = dl.TestFunction(self.Vh[STATE])
        up = dl.Function(self.Vh[STATE], x[STATE])
        a = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        res_varf = self._forwardvarf(self, a, up, test)
        Cvarf = dl.derivative(res_varf, a, trial)
        C = dl.assemble(Cvarf)
        return C
                
    def assembleWau(self, x):
        """
        Assemble the derivative of the parameter equation with respect to the state
        """
        trial = dl.TrialFunction(self.Vh[STATE])
        test  = dl.TestFunction(self.Vh[PARAMETER])
        vq = dl.Function(self.Vh[ADJOINT], x[ADJOINT])
        a = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        up = dl.Function(self.Vh[STATE], x[STATE])
        val = self._forwardvarf(self, a, up, vq)
        form = dl.derivative(val, a, test)
        varf = dl.derivative(form, up, trial)
        Wau = dl.assemble(varf)
        return Wau
    
    def assembleRaa(self, x):
        """
        Assemble the derivative of the parameter equation with respect to the parameter (Newton method)
        """
        trial = dl.TrialFunction(self.Vh[PARAMETER])
        test  = dl.TestFunction(self.Vh[PARAMETER])
        vq = dl.Function(self.Vh[ADJOINT], x[ADJOINT])
        a = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        up = dl.Function(self.Vh[STATE], x[STATE])
        val = self._forwardvarf(a, up, vq)
        form = dl.derivative(val,a,test)
        varf = dl.derivative(form, a, trial)
        return dl.assemble(varf)

            
    def cost(self, x):
        """
        Given the list x = [u,a,p] which describes the state, parameter, and
        adjoint variable compute the cost functional as the sum of 
        the misfit functional and the regularization functional.
        
        Return the list [cost functional, regularization functional, misfit functional]
        
        Note: p is not needed to compute the cost functional
        """        
        assert x[STATE] != None
                       
        diff = self.B*x[STATE]
        diff -= self.u_o
        misfit = (.5/self.noise_variance) * diff.inner(diff)
        
        Rdiff_x = dl.Vector()
        self.prior.init_vector(Rdiff_x,0)
        diff_x = x[PARAMETER] - self.prior.mean
        self.prior.R.mult(diff_x, Rdiff_x)
        reg = .5 * diff_x.inner(Rdiff_x)
        
        c = misfit + reg
        
        return c, reg, misfit
    
    def solveFwd(self, out, x, tol=1e-9):
        a = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        up = dl.Function(self.Vh[STATE], x[STATE])
        vq = dl.TestFunction(self.Vh[STATE])
        vh,qh = dl.split(vq)
        form = self._forwardvarf(a, up, vq)
        rhs  = self.f_vel*vh*dl.dx + dl.Constant(0.)*qh*dl.dx
        sol = dl.Function(self.Vh[STATE], out)
        dl.solve(form==rhs, sol, solver_parameters={"newton_solver": {"relative_tolerance":tol, "maximum_iterations":100}}))

    
    def solveAdj(self, out, x, tol=1e-9):
        """
        Solve the adjoint problem.
        """
        At = self.assembleA(x, assemble_adjoint = True)
        At.init_vector(out, 1)
        Bu = -(self.B*x[STATE])
        Bu += self.u_o
        rhs = dl.Vector()
        self.B.init_vector(rhs, 1)
        self.B.transpmult(Bu,rhs)
        rhs *= 1.0/self.noise_variance
        
        solver = dl.PETScLUSolver(At)
        nit = solver.solve(out,rhs)
    
    def evalGradientParameter(self,x, mg):
        """
        Evaluate the gradient for the variation parameter equation at the point x=[u,a,p].
        Parameters:
        - x = [u,a,p] the point at which to evaluate the gradient.
        - mg the variational gradient (g, atest) being atest a test function in the parameter space
          (Output parameter)
        
        Returns the norm of the gradient in the correct inner product g_norm = sqrt(g,g)
        """ 
        C = self.assembleC(x)

        self.prior.init_vector(mg,0)
        C.transpmult(x[ADJOINT], mg)
        Rdx = dl.Vector()
        self.prior.init_vector(Rdx,0)
        dx = x[PARAMETER] - self.prior.mean
        self.prior.R.mult(dx, Rdx)   
        mg.axpy(1., Rdx)
        
        g = dl.Vector()
        self.prior.init_vector(g,1)
        
        self.prior.Msolver.solve(g, mg)
        g_norm = dl.sqrt( g.inner(mg) )
        
        return g_norm
        
    
    def setPointForHessianEvaluations(self, x):  
        """
        Specify the point x = [u,a,p] at which the Hessian operator (or the Gauss-Newton approximation)
        need to be evaluated.
        """      
        self.A  = self.assembleA(x)
        self.At = self.assembleA(x, assemble_adjoint=True )
        self.C  = self.assembleC(x)
        self.Wau = self.assembleWau(x)
        self.Raa = self.assembleRaa(x)

        
    def solveFwdIncremental(self, sol, rhs, tol):
        """
        Solve the incremental forward problem for a given rhs
        """
        solver = dl.PETScLUSolver(self.A)
        self.A.init_vector(sol,1)
        nit = solver.solve(sol,rhs)
#        print "FwdInc", (self.A*sol-rhs).norm("l2")/rhs.norm("l2"), nit
        
    def solveAdjIncremental(self, sol, rhs, tol):
        """
        Solve the incremental adjoint problem for a given rhs
        """
        solver = dl.PETScLUSolver(self.At)
        self.At.init_vector(sol,1)
        nit = solver.solve(sol, rhs)
#        print "AdjInc", (self.At*sol-rhs).norm("l2")/rhs.norm("l2"), nit
    
    def applyC(self, da, out):
        self.C.mult(da,out)
    
    def applyCt(self, dp, out):
        self.C.transpmult(dp,out)
    
    def applyWuu(self, du, out):
        help = dl.Vector()
        self.B.init_vector(help, 0)
        self.B.mult(du, help)
        self.B.transpmult(help, out)
        out *= 1./self.noise_variance
    
    def applyWua(self, da, out):
        self.Wau.transpmult(da,out)

    
    def applyWau(self, du, out):
        self.Wau.mult(du, out)
    
    def applyR(self, da, out):
        self.prior.R.mult(da, out)
        
    def Rsolver(self):        
        return self.prior.Rsolver
    
    def applyRaa(self, da, out):
        self.Raa.mult(da, out)
            
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
    ntargets = 300
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
        pen = 1e1
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
    rel_noise = 0.01
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
    k = 50
    p = 20
    print "Double Pass Algorithm. Requested eigenvectors: {0}; Oversampling {1}.".format(k,p)
    Omega = np.random.randn(x[PARAMETER].array().shape[0], k+p)
    #d, U = singlePassG(Hmisfit, model.R, model.Rsolver, Omega, k, check_Bortho=True, check_Aortho=True, check_residual=True)
    d, U = doublePassG(Hmisfit, prior.R, prior.Rsolver, Omega, k, check_Bortho=False, check_Aortho=False, check_residual=False)
    posterior = GaussianLRPosterior(prior, d, U)
    posterior.mean = x[PARAMETER]
    
    post_tr, prior_tr, corr_tr = posterior.trace(method="Estimator", tol=1e-1, min_iter=20, max_iter=100)
    print "Posterior trace {0:5e}; Prior trace {1:5e}; Correction trace {2:5e}".format(post_tr, prior_tr, corr_tr)
    post_pw_variance, pr_pw_variance, corr_pw_variance = posterior.pointwise_variance("Exact")

    print sep, "Save State, Parameter, Adjoint, and observation in paraview", sep
    xxname = ["State", "Parameter", "Adjoint"]
    xx = [dl.Function(Vh[i], x[i], name=xxname[i]) for i in range(len(Vh))]
    dl.File("results/poisson_state.pvd") << xx[STATE]
    dl.File("results/poisson_state_true.pvd") << dl.Function(Vh[STATE], utrue, name = xxname[STATE])
    dl.File("results/poisson_parameter.pvd") << xx[PARAMETER]
    dl.File("results/poisson_parameter_true.pvd") << dl.Function(Vh[PARAMETER], atrue, name = xxname[PARAMETER])
    dl.File("results/poisson_parameter_prmean.pvd") << dl.Function(Vh[PARAMETER], prior.mean, name = xxname[PARAMETER])
    dl.File("results/poisson_adjoint.pvd") << xx[ADJOINT]
    
    vtrue = compute_velocity(mesh, Vh, atrue, utrue)
    dl.File("results/poisson_vel_true.pvd") << vtrue
    v_map = compute_velocity(mesh, Vh, x[PARAMETER], x[STATE])
    dl.File("results/poisson_vel.pvd") << v_map
    
    exportPointwiseObservation(targets, model.u_o, "results/poisson_observation.vtp")
    
    fid = dl.File("results/pointwise_variance.pvd")
    fid << dl.Function(Vh[PARAMETER], post_pw_variance, name="Posterior")
    fid << dl.Function(Vh[PARAMETER], pr_pw_variance, name="Prior")
    fid << dl.Function(Vh[PARAMETER], corr_pw_variance, name="Correction")
    
    
    print sep, "Generate samples from Prior and Posterior\n","Export generalized Eigenpairs", sep
    fid_prior = dl.File("samples/sample_prior.pvd")
    fid_post  = dl.File("samples/sample_post.pvd")
    nsamples = 500
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
        
    #Save eigenvalues for printing:
    posterior.exportU(Vh[PARAMETER], "hmisfit/evect.pvd")
    np.savetxt("hmisfit/eigevalues.dat", d)
    
    print sep, "Visualize results", sep
    dl.plot(xx[STATE], title = xxname[STATE])
    dl.plot(dl.exp(xx[PARAMETER]), title = xxname[PARAMETER])
    dl.plot(xx[ADJOINT], title = xxname[ADJOINT])
    
    plt.plot(range(0,k), d, 'b*', range(0,k), np.ones(k), '-r')
    plt.yscale('log')
        
    plt.show()    
    dl.interactive()
    
