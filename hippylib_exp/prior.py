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
import math

import sys
sys.path.append( "../")
from hippylib import *
from hippylib.prior import *

class ConstrainedBiLaplacianPrior(_Prior):
    """
    This class implement a Prior model with covariance matrix
    C = (\delta I + \gamma \div Theta \grad + \sum_i dirac(x - x_i) ) ^ {-2}.
    
    Here x_i (i=1,...,n) are points were we assume to know exactly the value
    of the parameter (i.e. m(x_i) = m_true( x_i) for i=1,...,n).
    
    The magnitude of \delta\gamma governs the variance of the samples, while
    the ratio \frac{\gamma}{\delta} governs the correlation lenght.
    
    Here Theta is a s.p.d tensor that models anisotropy in the covariance kernel.
    
    The prior mean is computed by solving the least square problem
    min || m ||^2_R, s.t. m(x_i) = m_true(x_i) for (i=1,...,n).
    """
    
    def __init__(self, Vh, gamma, delta, locations, m_true, Theta = None, pen = 1e4, rel_tol=1e-12, max_iter=1000):
        """
        Construct the Prior model.
        Input:
        - Vh:              the finite element space for the parameter
        - gamma and delta: the coefficient in the PDE
        - locations:       the points x_i at which we assume to know the
                           true value of the parameter
        - m_true:          the true model
        - Theta:           the s.p.d. tensor for anisotropic diffusion of the pde
        - pen:             a discretization parameter for the dirac delta
        """
        assert delta != 0. or pen != 0, "Intrinsic Gaussian Prior are not supported"
        self.Vh = Vh
        
        trial = dl.TrialFunction(Vh)
        test  = dl.TestFunction(Vh)
        
        if Theta == None:
            varfL = dl.inner(dl.nabla_grad(trial), dl.nabla_grad(test))*dl.dx
        else:
            varfL = dl.inner(Theta*dl.grad(trial), dl.grad(test))*dl.dx
        varfM = dl.inner(trial,test)*dl.dx
        
        self.M = dl.assemble(varfM)
        self.Msolver = dl.PETScKrylovSolver("cg", "jacobi")
        self.Msolver.set_operator(self.M)
        self.Msolver.parameters["maximum_iterations"] = max_iter
        self.Msolver.parameters["relative_tolerance"] = rel_tol
        self.Msolver.parameters["error_on_nonconvergence"] = True
        self.Msolver.parameters["nonzero_initial_guess"] = False
        
        Po = assemblePointwiseObservation(Vh, locations)
                
        self.A = dl.assemble(gamma*varfL+delta*varfM)
        PotPo = MatAtB(Po,Po)
        self.A.axpy(pen, PotPo, False);      
        self.Asolver = dl.PETScKrylovSolver("cg", amg_method())
        self.Asolver.set_operator(self.A)
        self.Asolver.parameters["maximum_iterations"] = max_iter
        self.Asolver.parameters["relative_tolerance"] = rel_tol
        self.Asolver.parameters["error_on_nonconvergence"] = True
        self.Asolver.parameters["nonzero_initial_guess"] = False
        
        Qh = dl.FunctionSpace(Vh.mesh(), 'Quadrature', 2*Vh._FunctionSpace___degree)
        ph = dl.TrialFunction(Qh)
        qh = dl.TestFunction(Qh)
        Mqh = dl.assemble(ph*qh*dl.dx)
        ones = dl.interpolate(dl.Constant(1.), Qh).vector()
        dMqh = Mqh*ones
        Mqh.zero()
        dMqh.set_local( ones.array() / np.sqrt(dMqh.array() ) )
        Mqh.set_diagonal(dMqh)
        MixedM = dl.assemble(ph*test*dl.dx)
        self.sqrtM = MatMatMult(MixedM, Mqh)
             
        self.R = _BilaplacianR(self.A, self.Msolver)      
        self.Rsolver = _BilaplacianRsolver(self.Asolver, self.M)
        
        rhs = dl.Vector()
        self.mean = dl.Vector()
        self.init_vector(rhs, 0)
        self.init_vector(self.mean, 0)
        
        PotPo.mult(m_true, rhs)
        rhs *= pen  
        self.Asolver.solve(self.mean, rhs)
        
     
    def init_vector(self,x,dim):
        """
        Inizialize a vector x to be compatible with the range/domain of R.
        If dim == "noise" inizialize x to be compatible with the size of
        white noise used for sampling.
        """
        if dim == "noise":
            self.sqrtM.init_vector(x, 1)
        else:
            self.A.init_vector(x,dim)   
        
    def sample(self, noise, s, add_mean=True):
        """
        Given a noise ~ N(0, I) compute a sample s from the prior.
        If add_mean=True add the prior mean value to s.
        """
        rhs = self.sqrtM*noise
        self.Asolver.solve(s, rhs)
        
        if add_mean:
            s.axpy(1., self.mean)


class LaplaceBeltramiPrior(_Prior):
    def __init__(self, Vh, gamma, delta, ds, mean=None, rel_tol=1e-12, max_iter=100):
        
        assert delta != 0., "Intrinsic Gaussian Prior are not supported"
        self.Vh = Vh
        trial = dl.TrialFunction(Vh)
        test  = dl.TestFunction(Vh)
        
        n = dl.FacetNormal(Vh.mesh()) 
        def dir_grad(uh,n):
            return dl.grad(uh) - dl.dot(dl.outer(n,n), dl.grad(uh))
        
        varfL = gamma*dl.inner(dir_grad(trial,n), dir_grad(test,n))*ds
        varfM = delta*dl.inner(trial,test)*ds
        
        self.M = dl.assemble(varfM)
        self.R = dl.assemble(varfL + varfM)
            
        self.M_1 = dl.assemble(varfM, keep_diagonal=True)
        self.R_1 = dl.assemble(varfL + varfM, keep_diagonal=True)
        self.M_1.ident_zeros()
        self.R_1.ident_zeros()
            
        self.Rsolver = dl.PETScKrylovSolver("cg", amg_method())
        self.Rsolver.set_operator(self.R_1)
        self.Rsolver.parameters["maximum_iterations"] = max_iter
        self.Rsolver.parameters["relative_tolerance"] = rel_tol
        self.Rsolver.parameters["error_on_nonconvergence"] = True
        self.Rsolver.parameters["nonzero_initial_guess"] = False
        
        self.Msolver = dl.PETScKrylovSolver("cg", "jacobi")
        self.Msolver.set_operator(self.M_1)
        self.Msolver.parameters["maximum_iterations"] = max_iter
        self.Msolver.parameters["relative_tolerance"] = rel_tol
        self.Msolver.parameters["error_on_nonconvergence"] = True
        self.Msolver.parameters["nonzero_initial_guess"] = False
        
        self.sampler = CGSampler()
        self.sampler.set_operator(self.R_1)
        
        self.mean = mean
        
        if self.mean is None:
            self.mean = dl.Vector()
            self.init_vector(self.mean, 0)
    
    def init_vector(self,x,dim):
        """
        Inizialize a vector x to be compatible with the range/domain of R.
        If dim == "noise" inizialize x to be compatible with the size of
        white noise used for sampling.
        """
        if dim == "noise":
            self.R.init_vector(x, 1)
        else:
            self.R.init_vector(x,dim)
            
    def sample(self, noise, s, add_mean=True):
        """
        Given a noise ~ N(0, I) compute a sample s from the prior.
        If add_mean=True add the prior mean value to s.
        """
        self.sampler.sample(noise.array(), s)
        
        if add_mean:
            s.axpy(1., self.mean)
