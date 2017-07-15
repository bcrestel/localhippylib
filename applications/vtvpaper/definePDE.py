"""
Define PDEs for both models
"""

from hippylib import PDEVariationalProblem
import dolfin as dl

def u_boundary(x, on_boundary):
    return on_boundary

def pdes(Vh, STATE, amg_method):
    f = dl.Constant(1.0)
    u_bdr = dl.Constant(0.0)
    u_bdr0 = dl.Constant(0.0)
    bc = dl.DirichletBC(Vh[STATE], u_bdr, u_boundary)
    bc0 = dl.DirichletBC(Vh[STATE], u_bdr0, u_boundary)
    
    def pde_varf(u,a,p):
        return dl.exp(a)*dl.inner(dl.nabla_grad(u), dl.nabla_grad(p))*dl.dx - f*p*dl.dx

    pde1 = PDEVariationalProblem(Vh, pde_varf, bc, bc0, is_fwd_linear=True)
    pde2 = PDEVariationalProblem(Vh, pde_varf, bc, bc0, is_fwd_linear=True)
    for pde in [pde1, pde2]:
        pde.solver = dl.PETScKrylovSolver("cg", amg_method())
        pde.solver.parameters["relative_tolerance"] = 1e-15
        pde.solver.parameters["absolute_tolerance"] = 1e-20
        pde.solver_fwd_inc = dl.PETScKrylovSolver("cg", amg_method())
        pde.solver_fwd_inc.parameters = pde.solver.parameters
        pde.solver_adj_inc = dl.PETScKrylovSolver("cg", amg_method())
        pde.solver_adj_inc.parameters = pde.solver.parameters

    return pde1, pde2
