# Copyright (c) 2016, The University of Texas at Austin & University of
# California, Merced.
#
# All Rights reserved.
# See file COPYRIGHT for details.
#
# This file is part of the hIPPYlib library. For more information and source code
# availability see https://hippylib.github.io.
#
# hIPPYlib is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.

import dolfin as dl
import sys
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

def prior_mean(Vh):
    cpp = \
    """
    class BetaPrior : public Expression
    {
    public:

    BetaPrior() :
    Expression()
    {
    }

void eval(Array<double>& values, const Array<double>& x) const
  {
  double val = 1000.;
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
    boundary_parts = dl.MeshFunction("size_t", mesh, "arolla_facet_region.xml")
    Vh2 = dl.VectorFunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
    Vh0 = dl.FunctionSpace(mesh, 'DG', 0)
    Vh_state =Vh2*Vh0
    Vh = [Vh_state, Vh1, Vh_state]
    print "Number of dofs: Velocity={0}, Pressure={1}".format(Vh2.dim(), Vh0.dim())
    
    ds = dl.Measure("ds")[boundary_parts]
    
    density = dl.Constant(910.) #kg/m^3
    gravity   = dl.Expression(("0.","-9.81")) #m/s^2
    GlenExp   = dl.Constant(3)
    RateFactor = dl.Constant(1e-16) #Pa^-n years^-1
    eps = dl.Constant(1e-5)
    pen = dl.Constant(1.e10)
    
    problem = IceProblem(Vh, ds(1), RateFactor, GlenExp, density, gravity, eps,pen)
    
    up_true = problem.generate_state()
    a_true = true_model(Vh1)
    
    problem.solveFwd(up_true, [up_true,a_true], mytol=None)
    
    up_true_fun = dl.Function(Vh[STATE], up_true)
    u_true_fun, p_true_fun = up_true_fun.split(deepcopy=True)
    MAX = u_true_fun.vector().norm("linf")
    print "MAX:", MAX
    
    gamma = 1.
    delta = 1.
    prior = LaplaceBeltramiPrior(Vh[PARAMETER], gamma, delta, ds(1), mean=prior_mean(Vh[PARAMETER]), rel_tol=1e-12, max_iter=100)
        
    u_trial,p_trial = dl.TrialFunctions(Vh[STATE])
    u_test,p_test = dl.TestFunctions(Vh[STATE])
    
    form = dl.inner(u_trial, u_test)*ds(2)
    misfit = ContinuousStateObservation(Vh[STATE], ds(2), None, form)
    
    misfit.d.set_local( up_true.array() )
    
    rel_noise = 0.1
    noise = problem.generate_state()
    noise.set_local(np.random.randn(Vh[STATE].dim()))
    misfit.d.axpy(MAX*rel_noise, noise)
    misfit.noise_variance = MAX*rel_noise*MAX*rel_noise
    
    model = Model(problem, prior,misfit)
    modelVerify(model,a_true, 1e-6)
    
    a0 = prior.mean.copy()
    solver = ReducedSpaceNewtonCG(model)
    solver.parameters["rel_tolerance"] = 1e-9
    solver.parameters["abs_tolerance"] = 1e-12
    solver.parameters["max_iter"]      = 25
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["c_armijo"] = 1e-4
    solver.parameters["GN_iter"] = 5
    
    x = solver.solve(a0)
    
    if solver.converged:
        print "\nConverged in ", solver.it, " iterations."
    else:
        print "\nNot Converged"

    print "Termination reason: ", solver.termination_reasons[solver.reason]
    print "Final gradient norm: ", solver.final_grad_norm
    print "Final cost: ", solver.final_cost

    
    plt.show()