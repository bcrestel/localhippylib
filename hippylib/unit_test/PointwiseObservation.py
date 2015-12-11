from dolfin import *
import sys
sys.path.append( "../../" )
from hippylib import *
import numpy as np

if __name__ == "__main__":
    set_log_active(False)
    ndim = 2
    nx = 64
    ny = 64
    mesh = UnitSquareMesh(nx, ny)
    
    ntargets = 3
    np.random.seed(seed=1)
    targets = np.random.uniform(0.1,0.9, [ntargets, ndim] )
    
    Vh = FunctionSpace(mesh, 'Lagrange', 1)
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
    for i in range(o.array().shape[0]):
        assert np.abs( o.array()[i] - targets[i,0] ) < 1e-10
    
    
    Vh2 = VectorFunctionSpace(mesh, 'Lagrange', 1)
    B2 = assemblePointwiseObservation(Vh2,targets)
    u2 = Vector()
    B2.init_vector(u2,1)
    
    o2 = Vector()
    B2.init_vector(o2,0)
    
    u2h = Function(Vh2, u2)
    u2h.interpolate(Expression(("x[0]", "x[1]") ))
        
    B2.mult(u2,o2)
    
    print o2.array()
    o2.set_local(o2.array())
    for i in range(targets.shape[0]):
        assert np.abs( o2.array()[2*i] - targets[i,0] ) < 1e-10
        assert np.abs( o2.array()[2*i+1] - targets[i,1] ) < 1e-10
        
    Xh = Vh2*Vh
    B3 = assemblePointwiseObservation(Xh,targets)
    up = Vector()
    B3.init_vector(up,1)
    
    o_up = Vector()
    B3.init_vector(o_up,0)
    
    uph = Function(Xh, up)
    uph.interpolate(Expression(("x[0]", "x[1]", "2.*x[0]+3.*x[1]+10." )))
        
    B3.mult(up,o_up)
    
    print o_up.array()
    for i in range(targets.shape[0]):
        assert np.abs( o_up.array()[3*i] - targets[i,0] ) < 1e-10
        assert np.abs( o_up.array()[3*i+1] - targets[i,1] ) < 1e-10
        assert np.abs( o_up.array()[3*i+2] - (2.*targets[i,0] + 3.*targets[i,1] + 10.) ) < 1e-10
    