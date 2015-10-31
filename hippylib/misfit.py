import dolfin as dl
from assemblePointwiseObservation import assemblePointwiseObservation
from variables import STATE
class Misfit:
    def cost(self,x):
        """Given x evaluate the cost functional"""
        
    def adj_rhs(self,x,rhs):
        """Evaluate the RHS for the adjoint problem"""
    
    def setLinearizationPoint(self,x):
        """Set the point for linearization. """
        
    def apply_ij(self,i,j, dir, out):
        """ apply \delta_ij cost in direction dir. """
        
class PointwiseStateObservation(Misfit):
    def __init__(self, Vh, obs_points):
        self.B = assemblePointwiseObservation(Vh, obs_points)
        self.d = dl.Vector()
        self.B.init_vector(self.d, 0)
        self.Bu = dl.Vector()
        self.B.init_vector(self.Bu, 0)
        self.noise_variance = 0.
        
    def cost(self,x):
        self.B.mult(x[STATE], self.Bu)
        self.Bu.axpy(-1., self.d)
        return (.5/self.noise_variance)*self.Bu.inner(self.Bu)
    
    def adj_rhs(self, x, out):
        self.B.mult(x[STATE], self.Bu)
        self.Bu.axpy(-1., self.d)
        self.B.transpmult(self.Bu, out)
        out *= (-1./self.noise_variance)
    
    def setLinearizationPoint(self,x):
        return
       
    def apply_ij(self,i,j,dir,out):
        if i == STATE and j == STATE:
            self.B.mult(dir, self.Bu)
            self.B.transpmult(self.Bu, out)
            out *= (1./self.noise_variance)
        else:
            out.zero()
            
class DistributedStateObservation(Misfit):
    def __init__(self, Vh, dX, bc, form = None):
        if form is None:
            u, v = dl.TrialFunction(Vh), dl.TestFunction(Vh)
            self.W = dl.assemble(dl.inner(u,v)*dX)
        else:
            self.W = dl.assemble( form )
                
        dummy = dl.Vector()
        if bc is not None:
            self.W.init_vector(dummy,0)
            self.bc.zero_columns(self.W, dummy)
            self.bc.zero(self.W)
        self.d = dl.Vector()
        self.W.init_vector(self.d,1)
        self.noise_variance = 0
        
    def cost(self,x):
        r = self.d.copy()
        r.axpy(-1., x[STATE])
        Wr = dl.Vector()
        self.W.init_vector(Wr,0)
        self.W.mult(r,Wr)
        return r.inner(Wr)/(2.*self.noise_variance)
    
    def adj_rhs(self, x, out):
        r = self.d.copy()
        r.axpy(-1., x[STATE])
        self.W.mult(r, out)
        out *= (1./self.noise_variance)
       
    def apply_ij(self,i,j,dir,out):
        if i == STATE and j == STATE:
            self.W.mult(dir, out)
            out *= (1./self.noise_variance)
        else:
            out.zero() 