from variables import STATE, PARAMETER, ADJOINT
from reducedHessian import ReducedHessian

import numpy as np
import matplotlib.pyplot as plt

def modelVerify(model,a0, innerTol):
    """
    Verify the reduced Gradient and the Hessian of a model.
    It will produce two loglog plots of the finite difference checks
    for the gradient and for the Hessian.
    It will also check for symmetry of the Hessian.
    """
    
    h = model.generate_vector(PARAMETER)
    h.set_local(np.random.normal(0, 1, len( h.array() )) )
    
    x = model.generate_vector()
    x[PARAMETER] = a0
    model.solveFwd(x[STATE], x, innerTol)
    model.solveAdj(x[ADJOINT], x, innerTol)
    cx = model.cost(x)
    
    grad_x = model.generate_vector(PARAMETER)
    model.evalGradientParameter(x, grad_x)
    grad_xh = grad_x.inner( h )
    
    model.setPointForHessianEvaluations(x)
    H = ReducedHessian(model,innerTol)
    Hh = model.generate_vector(PARAMETER)
    H.mult(h, Hh)
    
    n_eps = 32
    eps = np.power(.5, np.arange(n_eps))
    err_grad = np.zeros(n_eps)
    err_H = np.zeros(n_eps)
    
    for i in range(n_eps):
        my_eps = eps[i]
        
        x_plus = model.generate_vector()
        x_plus[PARAMETER].set_local( a0.array() )
        x_plus[PARAMETER].axpy(my_eps, h)
        model.solveFwd(x_plus[STATE],   x_plus, innerTol)
        model.solveAdj(x_plus[ADJOINT], x_plus,innerTol)
        
        dc = model.cost(x_plus)[0] - cx[0]
        err_grad[i] = abs(dc/my_eps - grad_xh)
        
        #Check the Hessian
        grad_xplus = model.generate_vector(PARAMETER)
        model.evalGradientParameter(x_plus, grad_xplus)
        
        err  = grad_xplus - grad_x
        err *= 1./my_eps
        err -= Hh
        
        err_H[i] = err.norm('linf')
        
    plt.figure()
    plt.subplot(121)
    plt.loglog(eps, err_grad, "-ob", eps, eps*(err_grad[0]/eps[0]), "-.k")
    plt.title("FD Gradient Check")
    plt.subplot(122)
    plt.loglog(eps, err_H, "-ob", eps, eps*(err_H[0]/eps[0]), "-.k")
    plt.title("FD Hessian Check")
    
        
    xx = model.generate_vector(PARAMETER)
    xx.set_local( np.random.normal(0, 1, len( xx.array() )) )
    yy = model.generate_vector(PARAMETER)
    yy.set_local( np.random.normal(0, 1, len( yy.array() )) )
    
    ytHx = H.inner(yy,xx)
    xtHy = H.inner(xx,yy)
    rel_symm_error = 2*abs(ytHx - xtHy)/(ytHx + xtHy)
    print "(yy, H xx) - (xx, H yy) = ", rel_symm_error
    if(rel_symm_error > 1e-10):
        print "HESSIAN IS NOT SYMMETRIC!!"


