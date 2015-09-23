'''
Created on Sep 17, 2015

@author: uvilla
'''
import dolfin as dl
import math
from hippylib import *

class GaussianDistribution:
    def __init__(self, H, a0):
        self.H = H
        self.a0 = a0
  
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
    