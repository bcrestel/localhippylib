import dolfin as dl
import sys
from ufl.operators import nabla_grad
from ufl.objects import dx
from dolfin.functions.function import TrialFunction
sys.path.append( "../../" )
from hippylib import *
import numpy as np
import matplotlib.pyplot as plt

from nonLinearStokesFwd import IceProblem

def true_model(Vh):
    return dl.interpolate(dl.Constant(1.), Vh).vector()
            
if __name__ == "__main__":
    dl.set_log_active(False)
    sep = "\n"+"#"*80+"\n"
    print sep, "Set up the mesh and finite element spaces", sep
    ndim = 2
    nx = 64
    ny = 64
    mesh = dl.Mesh("simple.xml")
    boundary_parts = dl.MeshFunction("size_t", mesh, "simple_facet_region.xml")
    Vh2 = dl.VectorFunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
    Vh0 = dl.FunctionSpace(mesh, 'DG', 0)
    Vh = Vh2*Vh1
    print "Number of dofs: Velocity={0}, Pressure={1}".format(Vh2.dim(), Vh0.dim())
    
    ds = dl.Measure("ds")[boundary_parts]
    
    density = dl.Constant(1.) #kg/m^3
    gravity   = dl.Expression(("0.","-1.")) #m/s^2
    GlenExp   = dl.Constant(3)
    RateFactor = dl.Constant(1) #Pa^-n years^-1
    eps = dl.Constant(.1)
    
    problem = IceProblem([Vh, Vh1], ds(1), RateFactor, GlenExp, density, gravity, eps)
    
    up = dl.Function(Vh)    
    vq = dl.TestFunction(Vh)
    a = dl.Function(Vh1, true_model(Vh1))
    
    Fnewtonian = problem.FNewtonean(up, a, vq, dl.Constant(1.))
    dl.solve(Fnewtonian==0, up)
    print "Mass conservation = ", dl.assemble( dl.dot( up.sub(0), dl.FacetNormal(mesh))*ds(2))
    dl.File("velocity_lin.pvd") << up.sub(0)
    dl.File("pressure_lin.pvd") << up.sub(1)
    
    E0 = dl.assemble(problem.Energy(up, a))  
    F0 = dl.assemble( problem.F(up,a,vq) )
    delta_u_F = dl.derivative( problem.F(up,a,vq), up)
    J = dl.assemble(delta_u_F)
    
    ## NOTE this is a div-free direction, otherwise it will not work :)
    dir = dl.interpolate(dl.Expression( ("x[1]", "x[0]", "0.")), Vh).vector()
#    dl.DirichletBC(Vh.sub(0), dl.Expression(("0.", "0.")), boundary_parts, 1).apply(dir)

    F0dir = F0.inner(dir)
    delta_u_F_dir = J*dir
    
    print "Energy = {0:5e}, ||FO|| = {1:5e}, ||J0dir||_2 = {2:5e}".format( E0, F0.norm("l2"), delta_u_F_dir.norm("l2") )
    
    n_eps = 32
    my_eps = 1.*np.power(.5, np.arange(n_eps))
    err_grad =np.zeros(n_eps)
    err_H = np.zeros(n_eps)
    
    up_2 = dl.Function(Vh)
    
    for i in range(n_eps):
        up_2.assign(up)
        up_2.vector().axpy(my_eps[i], dir )
        
        Eplus = dl.assemble( problem.Energy(up_2, a))
        err_grad[i] = abs( (Eplus-E0)/my_eps[i] - F0dir )
        
        Fplus = dl.assemble( problem.F(up_2, a, vq) )
        Fplus.axpy(-1., F0)
        Fplus *= 1./my_eps[i]
        Fplus.axpy(-1., delta_u_F_dir)
        err_H[i] = Fplus.norm("l2")
        
    plt.figure()
    plt.subplot(121)
    plt.loglog(my_eps, err_grad, "-ob")
    plt.loglog(my_eps, .5*err_grad[0] * (my_eps/my_eps[0]), "--k" )
    plt.title("FD Check gradient")
    plt.subplot(122)
    plt.title("FD Check Hessian")
    plt.loglog(my_eps, err_H, "-ob")
    plt.loglog(my_eps, .5*err_H[0] * (my_eps/my_eps[0]), "--k" )
    plt.show()
    
    x = [up.vector(), a.vector()]
    problem.solveFwd(up.vector(), x)
    
    dl.File("velocity.pvd") << up.sub(0)
    dl.File("pressure.pvd") << up.sub(1)
