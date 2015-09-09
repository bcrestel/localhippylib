from variables import STATE, CONTROL, ADJOINT
from reducedHessian import ReducedHessian

import numpy as np

def modelVerify(model,a0, eps, tol):
    
    innerTol = 1e-6*tol
    x = model.generate_vector()
    x[CONTROL] = a0
    model.solveFwd(x[STATE], x, innerTol)
    model.solveAdj(x[ADJOINT], x, innerTol)
    cx = model.cost(x)[0]
    
    #check the gradient
    h = model.generate_vector(CONTROL)
    h.set_local(np.random.normal(0, 1, len( h.array() )) )
#    print h.norm("l2")
    
    x_plus = model.generate_vector()
    x_plus[CONTROL] = a0 + eps*h
    model.solveFwd(x_plus[STATE],   x_plus, innerTol)
    model.solveAdj(x_plus[ADJOINT], x_plus,innerTol)
    
    x_minus = model.generate_vector()
    x_minus[CONTROL] = a0 - eps*h
    model.solveFwd(x_minus[STATE],  x_minus, innerTol)
    model.solveAdj(x_minus[ADJOINT], x_minus, innerTol)
    
    dc = model.cost(x_plus)[0] - model.cost(x_minus)[0]
    grad_x = model.generate_vector(CONTROL)
    model.evalGradientControl(x, grad_x)
    grad_xh = grad_x.inner( h )
    
    rel_err_grad = abs(dc/eps - 2.*grad_xh)/cx
    print "abs( c(x_plus) - c(x_minus) - 2*eps*(grad_x, h) )/(eps*c(x)) = ", rel_err_grad
    if(rel_err_grad > tol):
        print "CHECK THE GRADIENT COMPUTATION!!"
        
    #Check the Hessian
    grad_xplus = model.generate_vector(CONTROL)
    grad_xminus = model.generate_vector(CONTROL)
    model.evalGradientControl(x_plus, grad_xplus)
    model.evalGradientControl(x_minus, grad_xminus)
    
    model.setPointForHessianEvaluations(x)
    H = ReducedHessian(model,innerTol)
    Hh = model.generate_vector(CONTROL)
    H.mult(h, Hh)
    err  = grad_xplus
    err -= grad_xminus
    err *= 1./(2.*eps)
    err -= Hh
    
    rel_err_H = err.norm('linf')/grad_x.norm('linf')
    
    
    print "|| (grad_xplus - grad_xmins)/2eps - H_x h ||_inf/(grad_x) = ", rel_err_H
    if(rel_err_H > tol):
        print "CHECK THE HESSIAN COMPUTATION"
    
    xx = model.generate_vector(CONTROL)
    xx.set_local( np.random.normal(0, 1, len( xx.array() )) )
    yy = model.generate_vector(CONTROL)
    yy.set_local( np.random.normal(0, 1, len( yy.array() )) )
    
    ytHx = H.inner(yy,xx)
    xtHy = H.inner(xx,yy)
    rel_symm_error = 2*abs(ytHx - xtHy)/(ytHx + xtHy)
    print "(yy, H xx) - (xx, H yy) = ", rel_symm_error
    if(rel_symm_error > tol):
        print "HESSIAN IS NOT SYMMETRIC!!"


