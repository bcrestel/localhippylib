from dolfin import Vector
import numpy as np

class TimeDependentVector:
    
    def __init__(self, times, tol=1e-10):
        self.nsteps = len(times)
        self.data = [];
        
        for i in range(self.nsteps):
            self.data.append( Vector() )
             
        self.times = times
        self.tol = tol
    
    def copy(self, other):
        self.nsteps = other.nsteps
        self.times = other.times
        self.tol = other.tol
        self.data = []
        
        for v in other.data:
            self.data.append( v.copy() )
        
    def initialize(self,M,dim):
        for d in self.data:
            M.init_vector(d,dim)
            d.zero()
            
    def randn_perturb(self,std_dev):
        for d in self.data:
            noise = std_dev * np.random.normal(0, 1, len(d.array()))
            d.set_local(d.array() + noise)
    
    def axpy(self, a, other):
        for i in range(self.nsteps):
            self.data[i].axpy(a,other.data[i])
        
    def zero(self):
        for d in self.data:
            d.zero()
            
    def store(self, u, t):
        i = 0
        while i < self.nsteps-1 and 2*t > self.times[i] + self.times[i+1]:
            i += 1
            
        assert abs(t - self.times[i]) < self.tol
        
        self.data[i].set_local( u.array() )
        
    def retrieve(self, u, t):
        i = 0
        while i < self.nsteps-1 and 2*t > self.times[i] + self.times[i+1]:
            i += 1
            
        assert abs(t - self.times[i]) < self.tol
        
        u.set_local( self.data[i].array() )
        
    def norm(self, time_norm, space_norm):
        assert time_norm == "linf"
        s_norm = 0
        for i in range(self.nsteps):
            tmp = self.data[i].norm(space_norm)
            if tmp > s_norm:
                s_norm = tmp
        
        return s_norm
        
                