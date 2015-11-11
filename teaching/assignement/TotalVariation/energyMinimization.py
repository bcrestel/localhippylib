# Initialization
from dolfin import *
import math
import numpy as np
import logging
import matplotlib.pyplot as plt
import nb
from unconstrainedMinimization import InexactNewtonCG

logging.getLogger('FFC').setLevel(logging.WARNING)
logging.getLogger('UFL').setLevel(logging.WARNING)
set_log_active(False)

# Load the mesh and define the finite element space
mesh = Mesh("circle.xml")
Vh = FunctionSpace(mesh, "Lagrange", 2)

u = Function(Vh)
u_hat = TestFunction(Vh)
u_tilde = TrialFunction(Vh)

# Define the energy functional
f = Expression("100*exp(-100*(x[0]*x[0] + x[1]*x[1]))")
beta = Constant(100)
k1 = Constant(.1)
k2 = Constant(1.)

E = Constant(.5)*beta*u*u*ds + \
    Constant(.5)*(k1 + k2*u*u)*inner(nabla_grad(u), nabla_grad(u))*dx - \
    f*u*dx

grad = beta*u*u_hat*ds + (k2*u*u_hat)*inner(nabla_grad(u), nabla_grad(u))*dx + \
       (k1 + k2*u*u)*inner(nabla_grad(u), nabla_grad(u_hat))*dx - f*u_hat*dx

H = beta*u_tilde*u_hat*ds + \
    k2*u_tilde*u_hat*inner(nabla_grad(u), nabla_grad(u))*dx + \
    Constant(2.)*(k2*u*u_hat)*inner(nabla_grad(u_tilde), nabla_grad(u))*dx + \
    Constant(2.)*k2*u_tilde*u*inner(nabla_grad(u), nabla_grad(u_hat))*dx + \
    (k1 + k2*u*u)*inner(nabla_grad(u_tilde), nabla_grad(u_hat))*dx
    
    
solver = InexactNewtonCG()
solver.parameters["rel_tolerance"] = 1e-9
solver.parameters["abs_tolerance"] = 1e-12
solver.parameters["gdu_tolerance"] = 1e-18
solver.parameters["max_iter"] = 1000
solver.parameters["c_armijo"] = 1e-5
solver.parameters["print_level"] = 1
solver.parameters["max_backtracking_iter"] = 10

solver.solve(E, u, grad, H)

plot(u)
interactive()
