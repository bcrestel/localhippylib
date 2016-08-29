import math
from variables import STATE, PARAMETER, ADJOINT
import numpy as np
import dolfin as dl

class SampleStruct:
    def __init__(self, model, derivative_info=0):
        self.derivative_info = derivative_info
        self.m = model.generate_vector(PARAMETER)
        self.cost = 0
        
        if self.derivative_info >= 1:
            self.g = model.generate_vector(PARAMETER)
        else:
            self.g = None
        
        
    def assign(self, other):
        assert self.derivative_info == other.derivative_info
        self.cost = other.cost
        self.m = other.m.copy()
        
        if self.derivative_info >= 1:
            self.g = other.g.copy()

class MALA:
    def __init__(self, model):
        self.model = model
        self.parameters = {}
        self.parameters["inner_rel_tolerance"]   = 1e-9
        self.parameters["print_level"]           = 0
        self.parameters["number_of_samples"]     = 2000
        self.parameters["burn_in"]               = 1000
        self.parameters["delta_t"]               = 0.25*1e-4
        self.parameters["print_progress"]        = 20
        
        self.u = self.model.generate_vector(STATE)
        self.p = self.model.generate_vector(ADJOINT)

    def run(self, x_map, prior):
        number_of_samples = self.parameters["number_of_samples"]
        burn_in = self.parameters["burn_in"]
        inner_tol = self.parameters["inner_rel_tolerance"]

        current = SampleStruct(self.model, derivative_info=1)
        proposed = SampleStruct(self.model, derivative_info=1)
        
        current.m.zero()
        current.m.axpy(1., x_map[PARAMETER])
        self.model.solveFwd(self.u, [self.u,current.m,self.p], inner_tol)
        current.cost = self.model.cost([self.u,current.m,self.p])[2]
        self.model.solveAdj(self.p, [self.u,current.m,self.p], inner_tol)
        self.model.evalGradientParameter([self.u,current.m,self.p], current.g, misfit_only=True)
        
        print "Burn {0} samples".format(burn_in)
        sample_count = 0
        naccept = 0
        n_check = burn_in // self.parameters["print_progress"]
        while (sample_count < burn_in):
            naccept +=self.sample(current, proposed, prior)
            sample_count += 1
            if sample_count % n_check == 0:
                print "{0:2.1f} % completed, Acceptance ratio {1:2.1f} %".format(float(sample_count)/float(burn_in)*100,
                                                                         float(naccept)/float(sample_count)*100 )
        
        print "Generate {0} samples".format(number_of_samples)
        sample_count = 0
        naccept = 0
        n_check = number_of_samples // self.parameters["print_progress"]
        while (sample_count < number_of_samples):
            naccept +=self.sample(current, proposed, prior)
            sample_count += 1
            if sample_count % n_check == 0:
                print "{0:2.1f} % completed, Acceptance ratio {1:2.1f} %".format(float(sample_count)/float(number_of_samples)*100,
                                                                         float(naccept)/float(sample_count)*100 )

        
        return naccept
        
        
    def sample(self, current, proposed, prior): 
        inner_tol = self.parameters["inner_rel_tolerance"]
        proposed.m = self.proposal(current, prior)
        self.model.solveFwd(self.u, [self.u, proposed.m, self.p], )
        self.model.solveAdj(self.p, [self.u,proposed.m,self.p], inner_tol)
        proposed.cost = self.model.cost([self.u,proposed.m,self.p])[2]
        self.model.evalGradientParameter([self.u, proposed.m, self.p], proposed.g, misfit_only=True)
        rho_mp = self.acceptance_ratio(current, proposed, prior)
        rho_pm = self.acceptance_ratio(proposed, current, prior)
        al = rho_mp - rho_pm
        if(al > math.log(np.random.rand())):
            current.assign(proposed)
            return 1
        else:
            return 0

    def proposal(self, current, prior):
        delta_t = self.parameters["delta_t"]
        gradient_term = dl.Vector()
        prior.Rsolver.solve(gradient_term, current.g)
        noise = dl.Vector()
        prior.init_vector(noise, "noise")
        noise_size = noise.array().shape[0]
        noise.set_local( np.random.randn( noise_size ) )
        w = dl.Vector()
        prior.init_vector(w, 0)
        prior.sample(noise,w)
        delta_tp2 = 2 + delta_t
        d_gam = (2-delta_t)/(2+delta_t) * current.m  - (2*delta_t)/(delta_tp2)*gradient_term + math.sqrt(8*delta_t)/delta_tp2 * w
        return d_gam

    def acceptance_ratio(self, origin, destination, prior):
        delta_t = self.parameters["delta_t"]
        m_m = destination.m - origin.m
        p_m = destination.m + origin.m
        priori = dl.Vector()
        prior.Rsolver.solve(priori, origin.g)
        temp = priori.inner(origin.g)
        rho_uv = origin.cost + 0.5*origin.g.inner(m_m) + \
                0.25*delta_t*origin.g.inner(p_m) + \
                0.25*delta_t*temp
        return rho_uv
        




