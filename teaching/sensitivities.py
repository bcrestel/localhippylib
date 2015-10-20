# add description
# Import dependencies
from dolfin import *
import numpy as np

# create mesh and define function spaces
nx = 32
ny = 32
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, 'Lagrange', 1)

# define function for state and adjoint
u = Function(V)
p = Function(V)

# define Trial and Test Functions
du_trial, p_trial, u_trial = TrialFunction(V), TrialFunction(V), TrialFunction(V)
du_test,  p_test,  u_test  = TestFunction(V), TestFunction(V), TestFunction(V)

# initialize input functions
f  = Constant("1.0")
u0 = Constant("0.0")
m  = Constant("1.0")
dm = Expression('1.0')

# set up dirichlet boundary conditions
def u0_boundary(x,on_boundary):
    return on_boundary
bc = DirichletBC(V, u0, u0_boundary)

# variational forms
var_state     = inner( m * nabla_grad(u_trial), nabla_grad(u_test)) * dx
var_rhs_state = f * u_test * dx
var_mass      = inner(u_trial, u_test) * dx

mat_state, rhs_state = assemble_system(var_state, var_rhs_state, bc)


mass = assemble(var_mass)

# solve state equation
u = Function(V)
solve(mat_state, u.vector(), rhs_state)

if mesh.num_cells() < 16: # print for small meshes only
    print mat_state.array()
    print mass.array()
    print rhs_state.array()

u.rename("state","ignore this")
File("state.pvd") << u

# solve sensitivity equations
var_sens      = inner( m * nabla_grad(du_trial), nabla_grad(du_test)) * dx
var_rhs_sens  = inner( -dm * nabla_grad(u), nabla_grad(du_test) ) * dx

mat_sens, rhs_sens   = assemble_system(var_sens, var_rhs_sens, bc)

du = Function(V,name="sens1")
solve(mat_sens, du.vector(), rhs_sens)
File("sens1.pvd") << du

## misfit
ud = Function(V,name="data")
ud.assign(u)
## perturb state solution and create synthetic measurements ud
## ud = u + ||u||/SNR * random.normal
noise_level = 0.05
MAX         = ud.vector().norm("linf")
noise       = Vector()
mat_sens.init_vector(noise,1)
noise.set_local( noise_level * MAX * np.random.normal(0, 1, len(ud.vector().array())) )
bc.apply(noise)
ud.vector().axpy(1., noise)
File("data.pvd") << ud

var_misfit = inner( (u - ud), u_test)*dx
vec_misfit = assemble(var_misfit)
grad1 = vec_misfit.inner(du.vector())
print grad1

# adjoint approach
var_adj = inner( m * nabla_grad(p_trial), nabla_grad(p_test)) * dx
rhs_adj = -var_misfit

p = Function(V,name="adj")
solve(var_adj==rhs_adj,p,bc)
File("adj.pvd") << p
grad_adj = assemble(inner (dm * nabla_grad(u), nabla_grad(p) )* dx)
print grad_adj
