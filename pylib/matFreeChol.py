from dolfin import Vector
import numpy as np

class MatFreeChol:
    """
    For the details regarding the A-orthogonalization algorithm see:
    
    Michele Benzi and Miroslav Tuma
    A robust incomplete factorization preconditioner for positive definite matrices
    Numer. Linear Algebra Appl. 2003; 10:385-400 (DOI: 10.1002/nla.320)
        
    This class should not be used as a solver O(n^3), but more as a sampler.
    """
    
    def __init__(self, A):
        self.A = A
        n = A.size(0)
        self.Z = np.eye(n, dtype=np.float64)
        self.Dinv = np.zeros(n, dtype=np.float64)
        
        zj, Azj = Vector(), Vector()
        self.A.init_vector(zj, 0)
        self.A.init_vector(Azj, 0)
        
        # this is O(n^3) but it can be optimized
        for j in range(0,n):
            zj.set_local( self.Z[:,j] )
            self.A.mult(zj, Azj)
            dinv = 1. /zj.inner(Azj)
            self.Dinv[j] = dinv
            for i in range(j+1,n):
                zi = self.Z[:,i]
                s = np.dot(zi.T, Azj.array())*dinv
                self.Z[:,i] = zi - s*zj.array()
                
    def init_vector(self,x, dim):
        self.A.init_vector(x)
        
    def solve(self, sol, rhs):
        Ztrhs = np.dot( self.Z.T, rhs.array() )
        DinvZtrhs = self.Dinv * Ztrhs
        sol.set_local( np.dot( self.Z, DinvZtrhs ) )
        
    def solve_sqrt(self, sol, rhs):
        Ztrhs = np.dot( self.Z.T, rhs.array() )
        sol.set_local(np.sqrt(self.Dinv) * Ztrhs)
        
    def solve_sqrt_t(self, sol, rhs):
        dinvrhs = np.sqrt(self.Dinv) * rhs.array()
        sol.set_local( np.dot(self.Z, dinvrhs))