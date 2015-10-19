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
    cpp = \
    """
    class BetaTrue : public Expression
    {
    public:

    BetaTrue() :
    Expression()
    {
    }

void eval(Array<double>& values, const Array<double>& x) const
  {
  double pi = 4.*atan(1.);
  double val = 0.;
  double x0 = x[0];
  if(x0 < 3750.)
     val = 1000.*(1. + sin(2*pi*x0/5000.));
  else if (x0 < 4000.)
     val = 1000.*(16 - x0/250.);
  else
     val = 1000.;
  values[0] = log(val);
  }  
};
    """
    return dl.interpolate(dl.Expression(cpp), Vh).vector()
            
if __name__ == "__main__":
    dl.set_log_active(False)
    sep = "\n"+"#"*80+"\n"
    print sep, "Set up the mesh and finite element spaces", sep
    ndim = 2
    nx = 64
    ny = 64
    mesh = dl.Mesh("arolla.xml")
    #boundary_parts = dl.FacetFunction("size_t", mesh)
    boundary_parts = dl.MeshFunction("size_t", mesh, "arolla_facet_region.xml")
#    dl.File("arolla_facet_region.xml") >> boundary_parts
    Vh2 = dl.VectorFunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
    Vh0 = dl.FunctionSpace(mesh, 'DG', 0)
    Vh = Vh2*Vh1
    print "Number of dofs: Velocity={0}, Pressure={1}".format(Vh2.dim(), Vh0.dim())
    
    ds = dl.Measure("ds")[boundary_parts]
    
    density = dl.Constant(910.) #kg/m^3
    gravity   = dl.Expression(("0.","-9.81")) #m/s^2
    GlenExp   = dl.Constant(3)
    RateFactor = dl.Constant(1e-16) #Pa^-n years^-1
    eps = dl.Constant(10.)
    
    problem = IceProblem(Vh, ds(1), RateFactor, GlenExp, density, gravity, eps)
    
    up = dl.Function(Vh)    
    vq = dl.TestFunction(Vh)
    a = dl.Function(Vh1, true_model(Vh1))
    
    Fnewtonian = problem.FNewtonean(up, a, vq, dl.Constant(1.))
    dl.solve(Fnewtonian==0, up)
       
    F0 = dl.assemble( problem.F(up,a,vq) )
    print "||FO|| = ", F0.norm("l2")
    delta_u_F = dl.derivative( problem.F(up,a,vq), up)
    J = dl.assemble(delta_u_F)
#    up_trial = dl.TrialFunction(Vh)    
#    J = dl.assemble( problem.F(up_trial, a, vq) )
    
    dir = dl.interpolate(dl.Expression( ("1e3*log(x[1]+1.)", "1e3*log(x[0]+1.)", "0.")), Vh).vector()
    dl.DirichletBC(Vh.sub(0), dl.Expression(("0.", "0.")), boundary_parts, 1).apply(dir)

    delta_u_F_dir = J*dir
    print "||Jdir||_2 = ", delta_u_F_dir.norm("l2")
    
    n_eps = 32
    my_eps = 1.*np.power(.5, np.arange(n_eps))
    err =np.zeros(n_eps)
    
    up_2 = dl.Function(Vh)
    for i in range(n_eps):
        up_2.assign(up)
        up_2.vector().axpy(my_eps[i], dir )
        Fplus = dl.assemble( problem.F(up_2, a, vq) )
        a1 =  Fplus.norm("l2")
        Fplus.axpy(-1., F0)
        a2 = Fplus.norm("l2")
        Fplus *= 1./my_eps[i]
        a3 = Fplus.norm("l2")
        Fplus.axpy(-1., delta_u_F_dir)
        print a1, a2, a3
        err[i] = Fplus.norm("l2")
        
    plt.figure()
    plt.loglog(my_eps, err, "-ob")
    plt.loglog(my_eps, .5*err[0] * (my_eps/my_eps[0]), "--k" )
    plt.show()
        
        
    dl.File("velocity.pvd")<< up.sub(0)
    dl.File("pressure.pvd") << up.sub(1)
