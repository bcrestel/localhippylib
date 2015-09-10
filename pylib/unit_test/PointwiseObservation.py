from dolfin import *
import sys
sys.path.append( "../../" )
from pylib import *
import numpy as np

if __name__ == "__main__":
    set_log_active(False)
    ndim = 2
    nx = 64
    ny = 64
    mesh = UnitSquareMesh(nx, ny)
    Vh = FunctionSpace(mesh, 'Lagrange', 1)
    
    ntargets = 3
    np.random.seed(seed=1)
    targets = np.random.uniform(0.1,0.9, [ntargets, ndim] )
    
    B = assemblePointwiseObservation(Vh,targets)
    
    u = Vector()
    B.init_vector(u,1)
    
    o = Vector()
    B.init_vector(o,0)
    
    uh = Function(Vh, u)
    uh.interpolate(Expression("x[0]") )
        
    B.mult(u,o)
    
    print targets
    print o.array()
    
    o.set_local(o.array())
    
    plot(uh)
    interactive()
    
    
    