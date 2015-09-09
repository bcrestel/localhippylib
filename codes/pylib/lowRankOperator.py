from dolfin import Vector
import numpy as np

class LowRankOperator:
    """ model the action of A = U D U^T.
    Where D is a diagonal matrix.
    
    """
    def __init__(self,d,U):
        self.d = d
        self.U = U
        
    def init_vector(self, x, dim):
        raise NameError("Not Implemented")
        
    def mult(self,x,y):
        Utx = np.dot( self.U.T, x.array() )
        dUtx = self.d*Utx
        y.set_local(np.dot(self.U, dUtx))
        
    def solve(self, sol, rhs):
        Utx = np.dot( self.U.T, rhs.array() )
        dinvUtx = Utx / self.d
        sol.set_local(np.dot(self.U, dinvUtx))
        
    def get_diagonal(self, diag):
        V = self.U * self.d
        diag.set_local(np.sum(V*self.U, 1))
        
    def trace(self,W=None):
        if W is None:
            diagUtU = np.sum(self.U*self.U,0)
            tr = np.sum(self.d*diagUtU)
        else:
            WU = np.zeros(self.U.shape, dtype=self.U.dtype)
            u, wu = Vector(), Vector()
            W.init_vector(u,1)
            W.init_vector(wu,0)
            for i in range(self.U.shape[1]):
                u.set_local(self.U[:,i])
                W.mult(u,wu)
                WU[:,i] = wu.array()
            diagWUtU = np.sum(WU*self.U,0)
            tr = np.sum(self.d*diagWUtU)
            
        return tr