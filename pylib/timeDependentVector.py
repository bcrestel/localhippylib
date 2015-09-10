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
        
    def initialize(self,M,dim):
        for d in self.data:
            M.init_vector(d,dim)
            d.zero()
    
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
        
                