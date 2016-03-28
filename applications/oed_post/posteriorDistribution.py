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
import math
from hippylib import *

class GaussianDistribution:
    def __init__(self, H, a0):
        self.H = H
        self.a0 = a0
        
    def generate_vector(self):
        out = dl.Vector()
        self.H.init_vector(out, 0)
        return out
  
    def __call__(self,a):
        r = a - self.a0
        c = self.H.inner(r,r)
        return math.exp(-.5*c)
    
    def gradient(self, a0, out):
        r = a0 - self.a0
        self.H.mult(r,out)
        c = .5*out.inner(r)
        out *= -math.exp(-c)
        
    def hessian_apply(self, a0, da, dout):
        r = a0 - self.a0
        grad = dl.Vector()
        self.H.init_vector(grad, 0)
        self.H.mult(r, grad)
        c = r.inner(grad)
        self.H.mult(da, dout)
        grad_da = grad.inner(da)
        dout.axpy(-grad_da, grad)
        dout *= -math.exp(-.5*c)
        
                
class PosteriorDistribution:
    def __init__(self, model, cmap):
        self.model = model
        self.x   = self.model.generate_vector()
        self.mg = self.model.generate_vector(PARAMETER)
        self.H = ReducedHessian(self.model, 1e-15)
        self.cmap = cmap
        
    def generate_vector(self):
        out = dl.Vector()
        self.H.init_vector(out, 0)
        return out
        
    def __call__(self, a):
        self.x[PARAMETER] = a.copy()
        self.model.solveFwd(self.x[STATE], self.x, 1e-15)
        [c, reg, misfit] = self.model.cost(self.x)
        return math.exp( -(c-self.cmap) )
    
    def gradient(self, a0, out):
        self.x[PARAMETER] = a0.copy()
        self.model.solveFwd(self.x[STATE], self.x, 1e-15)
        self.model.solveAdj(self.x[ADJOINT], self.x, 1e-15)
        [c, reg, misfit] = self.model.cost(self.x)
        self.model.evalGradientParameter(self.x, out)
        out *= -math.exp( -(c-self.cmap) )
        
    def hessian_apply(self, a0, da, dout):
        self.x[PARAMETER] = a0.copy()
        self.model.solveFwd(self.x[STATE], self.x, 1e-15)
        self.model.solveAdj(self.x[ADJOINT], self.x, 1e-15)
        [c, reg, misfit] = self.model.cost(self.x)
        self.model.evalGradientParameter(self.x, self.mg)
        self.model.setPointForHessianEvaluations(self.x)
        self.H.mult(da, dout)
        dout.axpy(-self.mg.inner(da), self.mg)
        dout *= -math.exp( -(c-self.cmap) )
    
def check_derivatives_pdf(pdf,start, dir):
                
    n = 32
    all_eps = 1e-1*np.power(2., -np.arange(0,n))
    err_grad = np.zeros(n)
    err_H = np.zeros(n)
    grad = pdf.generate_vector()
    pdf.gradient(start, grad)
    graddir = grad.inner(dir)
    c = pdf(start)
    for i in range(n):
        eps = all_eps[i]
        aplus = start.copy()
        aplus.axpy(eps, dir)        
        cplus = pdf(aplus)
        err_grad[i] =  np.abs( (cplus-c)/eps - graddir)

    
    Hdir = pdf.generate_vector()
    pdf.hessian_apply(start, dir, Hdir)
    grad_plus = pdf.generate_vector()
    for i in range(n):
        eps = all_eps[i]
        aplus = start.copy()
        aplus.axpy(eps, dir)        
        pdf.gradient(aplus, grad_plus)        
        err = grad_plus.copy()
        err.axpy(-1, grad)
        err *= 1./eps
        err.axpy(-1., Hdir)
        err_H[i] = err.norm("linf")
               
    plt.figure()
    plt.subplot(121)
    plt.loglog(all_eps, err_grad, "-ob", all_eps, all_eps, "--k")
    plt.ylabel("Error")
    plt.xlabel("eps")
    plt.title("Finite Difference check (Gradient)")
    plt.subplot(122)
    plt.loglog(all_eps, err_H,"-ob", all_eps, all_eps, "--k")
    plt.ylabel("Error")
    plt.xlabel("eps")
    plt.title("Finite Difference check (Hessian)")
    plt.show()


