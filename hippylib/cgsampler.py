from dolfin import Vector
import numpy as np
import math

class CGSampler:
    """ Implementation of the CG Sampler as in
        Albert Parker and Colin Fox
        Sampling Gaussian Distributions in Krylov Spaces with Conjugate Gradient
        SIAM J SCI COMPUT, Vol 34, No. 3 pp. B312-B334
    """

    def __init__(self):
        self.parameters = {}
        self.parameters["tolerance"] = 1e-4
        self.parameters["print_level"] = 0
        self.parameters["verbose"] = 0
        
        self.A = None
        self.converged = False
        self.iter = 0
        
        self.b = Vector()
        self.r = Vector()
        self.p = Vector()
        self.Ap = Vector()
                
    def set_operator(self, A):
        self.A = A
        self.A.init_vector(self.r,0)
        self.A.init_vector(self.p,0)
        self.A.init_vector(self.Ap,0)
        
        self.A.init_vector(self.b,0)
        self.b.set_local(np.random.randn(self.b.array().shape[0]))
                        
    def sample(self, noise, s):
        
        s.zero()
        
        self.iter = 0
        self.converged = False
        
        # r0 = b
        self.r.zero()
        self.r.axpy(1., self.b)
        
        #p0 = r0
        self.p.zero()
        self.p.axpy(1., self.r)
        
        self.A.mult(self.p, self.Ap)
        
        d = self.p.inner(self.Ap)
        
        tol2 = self.parameters["tolerance"]*self.parameters["tolerance"]
        
        rnorm2_old = self.r.inner(self.r)
        
        if self.parameters["verbose"] > 0:
            print "initial residual = ", math.sqrt(rnorm2_old)
        
        while (not self.converged) and (self.iter < noise.shape[0]):
            gamma = rnorm2_old/d
            s.axpy(noise[self.iter]/math.sqrt(d), self.p)
            self.r.axpy(-gamma, self.Ap)
            rnorm2 = self.r.inner(self.r)
            beta = rnorm2/rnorm2_old
            # p_new = r + beta p
            self.p *= beta
            self.p.axpy(1., self.r)
            self.A.mult(self.p, self.Ap)
            d = self.p.inner(self.Ap)
            rnorm2_old = rnorm2
            
            if rnorm2 < tol2:
                self.converged = True
            else:
                rnorm2_old = rnorm2
                self.iter = self.iter+1
         
        if self.parameters["verbose"] > 0:       
            print "Final residual {0} after {1} iterations".format( math.sqrt(rnorm2_old), self.iter)
            
        