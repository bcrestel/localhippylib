import math
from variables import STATE, PARAMETER, ADJOINT
import numpy as np
import dolfin as dl

class NullQoi(object):
    def __init__(self):
        pass
    def eval(self,x):
        return 0.

class NullTracer(object):
    def __init__(self):
        pass
    def append(self,current, q):
        pass
    
class QoiTracer(object):
    def __init__(self, n):
        self.qs = np.zeros(n)
        self.i = 0
        
    def append(self,current, q):
        self.qs[self.i] = q
        self.i+=1

class SampleStruct:
    def __init__(self, kernel):
        self.derivative_info = kernel.derivativeInfo()
        self.u = kernel.model.generate_vector(STATE)
        self.m = kernel.model.generate_vector(PARAMETER)
        self.cost = 0
        
        if self.derivative_info >= 1:
            self.p = kernel.model.generate_vector(STATE)
            self.g = kernel.model.generate_vector(PARAMETER)
        else:
            self.p = None
            self.g = None
        
    def assign(self, other):
        assert self.derivative_info == other.derivative_info
        self.cost = other.cost
        
        self.m = other.m.copy()
        self.u = other.u.copy()
        
        if self.derivative_info >= 1:
            self.g = other.g.copy()
            self.p = other.p.copy()


class MCMC(object):
    def __init__(self, kernel):
        self.kernel = kernel
        self.parameters = {}
        self.parameters["number_of_samples"]     = 2000
        self.parameters["burn_in"]               = 1000
        self.parameters["print_progress"]        = 20
        
        self.sum_q = 0.
        self.sum_q2 = 0.
        
    def run(self, m0, qoi=None, tracer = None):
        if qoi is None:
            qoi = NullQoi()
        if tracer is None:
            tracer = NullTracer()
        number_of_samples = self.parameters["number_of_samples"]
        burn_in = self.parameters["burn_in"]
        
        current = SampleStruct(self.kernel)
        proposed = SampleStruct(self.kernel)
        
        current.m.zero()
        current.m.axpy(1., m0)
        self.kernel.init_sample(current)
        
        print "Burn {0} samples".format(burn_in)
        sample_count = 0
        naccept = 0
        n_check = burn_in // self.parameters["print_progress"]
        while (sample_count < burn_in):
            naccept +=self.kernel.sample(current, proposed)
            sample_count += 1
            if sample_count % n_check == 0:
                print "{0:2.1f} % completed, Acceptance ratio {1:2.1f} %".format(float(sample_count)/float(burn_in)*100,
                                                                         float(naccept)/float(sample_count)*100 )
        
        print "Generate {0} samples".format(number_of_samples)
        sample_count = 0
        naccept = 0
        n_check = number_of_samples // self.parameters["print_progress"]
        while (sample_count < number_of_samples):
            naccept +=self.kernel.sample(current, proposed)
            q = qoi.eval([current.u, current.m])
            self.sum_q += q
            self.sum_q2 += q*q
            tracer.append(current, q)
            sample_count += 1
            if sample_count % n_check == 0:
                print "{0:2.1f} % completed, Acceptance ratio {1:2.1f} %".format(float(sample_count)/float(number_of_samples)*100,
                                                                         float(naccept)/float(sample_count)*100 )        
        return naccept


class MALAKernel:
    def __init__(self, model):
        self.model = model
        self.parameters = {}
        self.parameters["inner_rel_tolerance"]   = 1e-9
        self.parameters["delta_t"]               = 0.25*1e-4
        
    def derivativeInfo(self):
        return 1

    def init_sample(self, s):
        inner_tol = self.parameters["inner_rel_tolerance"]
        self.model.solveFwd(s.u, [s.u,s.m,s.p], inner_tol)
        s.cost = self.model.cost([s.u,s.m,s.p])[2]
        self.model.solveAdj(s.p, [s.u,s.m,s.p], inner_tol)
        self.model.evalGradientParameter([s.u,s.m,s.p], s.g, misfit_only=True) 
        
    def sample(self, current, proposed): 
        proposed.m = self.proposal(current)
        self.init_sample(proposed)
        rho_mp = self.acceptance_ratio(current, proposed)
        rho_pm = self.acceptance_ratio(proposed, current)
        al = rho_mp - rho_pm
        if(al > math.log(np.random.rand())):
            current.assign(proposed)
            return 1
        else:
            return 0

    def proposal(self, current):
        delta_t = self.parameters["delta_t"]
        gradient_term = dl.Vector()
        self.model.prior.Rsolver.solve(gradient_term, current.g)
        noise = dl.Vector()
        self.model.prior.init_vector(noise, "noise")
        noise_size = noise.array().shape[0]
        noise.set_local( np.random.randn( noise_size ) )
        w = dl.Vector()
        self.model.prior.init_vector(w, 0)
        self.model.prior.sample(noise,w)
        delta_tp2 = 2 + delta_t
        d_gam = (2-delta_t)/(2+delta_t) * current.m  - (2*delta_t)/(delta_tp2)*gradient_term + math.sqrt(8*delta_t)/delta_tp2 * w
        return d_gam

    def acceptance_ratio(self, origin, destination):
        delta_t = self.parameters["delta_t"]
        m_m = destination.m - origin.m
        p_m = destination.m + origin.m
        priori = dl.Vector()
        self.model.prior.Rsolver.solve(priori, origin.g)
        temp = priori.inner(origin.g)
        rho_uv = origin.cost + 0.5*origin.g.inner(m_m) + \
                0.25*delta_t*origin.g.inner(p_m) + \
                0.25*delta_t*temp
        return rho_uv
        

class pCNKernel:
    def __init__(self, model):
        self.model = model
        self.parameters = {}
        self.parameters["inner_rel_tolerance"]   = 1e-9
        self.parameters["s"]                     = 0.1

    def derivativeInfo(self):
        return 0

    def init_sample(self, current):
        inner_tol = self.parameters["inner_rel_tolerance"]
        self.model.solveFwd(current.u, [current.u,current.m,None], inner_tol)
        current.cost = self.model.cost([current.u,current.m,None])[2]
        
    def sample(self, current, proposed): 
        proposed.m = self.proposal(current)
        self.init_sample(proposed)
        al = -proposed.cost + current.cost
        if(al > math.log(np.random.rand())):
            current.assign(proposed)
            return 1
        else:
            return 0

    def proposal(self, current):
        #Generate sample from the prior
        noise = dl.Vector()
        self.model.prior.init_vector(noise, "noise")
        noise_size = noise.array().shape[0]
        noise.set_local( np.random.randn( noise_size ) )
        noise.apply("")
        w = dl.Vector()
        self.model.prior.init_vector(w, 0)
        self.model.prior.sample(noise,w, add_mean=False)
        # do pCN linear combination with current sample
        s = self.parameters["s"]
        w *= s
        w.axpy(1., self.model.prior.mean)
        w.axpy(np.sqrt(1. - s*s), current.m - self.model.prior.mean)
        
        return w





