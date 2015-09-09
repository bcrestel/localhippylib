#J. Freund, R. Stenberg. On weakly imposed boundary conditions for second order problems.
#Proceedings of the Ninth Int. Conf. Finite Elements in Fluids, Venice 1995.
#M. Morandi Cecchi et al., Eds. pp. 327-336

from dolfin import *
import sys
from dolfin.cpp.la import PETScKrylovSolver
sys.path.append( "../../" )
import numpy as np
from pylib import *

nx = 64
ny = 64

mesh = UnitSquareMesh(nx, ny)

class TopBoundary(SubDomain):
    def inside(self,x,on_boundary):
        y_in = x[1] > 1-DOLFIN_EPS and x[1] < 1 + DOLFIN_EPS
        return on_boundary and y_in
    
class LateralBoundary(SubDomain):
    def inside(self,x,on_boundary):
        x_left = x[0] > -DOLFIN_EPS and x[0] < DOLFIN_EPS
        x_right = x[0] > 1.-DOLFIN_EPS and x[0] < 1.+DOLFIN_EPS
        return on_boundary and ( x_left or x_right )

class LateralPeriodicBoundary(SubDomain):
    def inside(self,x,on_boundary):
        x_left = x[0] > -DOLFIN_EPS and x[0] < DOLFIN_EPS
        return on_boundary and x_left
    
    def map(self,x,y):
        y[0] = x[0] - 5.
        y[1] = x[1]       
    
class BottomBoundary(SubDomain):
    def inside(self,x,on_boundary):
        y_in = x[1] > -DOLFIN_EPS and x[1] < DOLFIN_EPS
        return on_boundary and y_in
    
boundary_markers = FacetFunction("size_t",mesh)

Gamma_B = BottomBoundary()
Gamma_B.mark(boundary_markers,1)

Gamma_T = TopBoundary()
Gamma_T.mark(boundary_markers,2)

Gamma_L = LateralPeriodicBoundary()

ds = Measure("ds",domain=mesh, subdomain_data=boundary_markers)

def coordinateTransform(xy):
    print xy.shape
    xynew = np.zeros(xy.shape, dtype=xy.dtype)
    xynew[:,0] = xy[:,0]*5
    #xynew[:,1] = (xy[:,1] -1. )*(.1*pow(xy[:,0],2)+1.) + 1.
    xynew[:,1] = xy[:,1]
    return xynew
    
mesh.coordinates()[:] = coordinateTransform( mesh.coordinates() )

Vh = VectorFunctionSpace(mesh, 'Lagrange', 2, constrained_domain=LateralPeriodicBoundary())
Wh = FunctionSpace(mesh, 'DG', 0, constrained_domain=LateralPeriodicBoundary())
Xh = MixedFunctionSpace([Vh,Wh])

beta = Expression("100.0")
eta = Expression("1.0")
rhog = Expression( ("g*sin(theta)", "-g*cos(theta)"), theta = 0.1*pi/180, g = 9.81 )

(uh, ph) = TrialFunctions(Xh)
(vh, qh) = TestFunctions(Xh)

def strain(uh):
    return sym(nabla_grad(uh))

def sigma_n(uh, ph, n):
    return -ph*n + Constant(2.)*eta*strain(uh)*n

def t(u,p,n): return dot(2.*eta*sym(grad(u)),n) - p*n
def tu(u,n): return dot(2.*eta*sym(grad(u)),n)
def T(u,n, beta): return beta*(uh - dot(outer(n,n),uh)) 

pen = Constant(100.)
n = FacetNormal(mesh)
h = CellSize(mesh)

stokesvarf = inner(Constant(2.)*eta*strain(uh), strain(vh))*dx - ph*nabla_div(vh)*dx - nabla_div(uh)*qh*dx
robinvarf  = inner( T(uh,n,beta), vh) * ds(1)
weak_bc = -dot(n, t(uh, ph, n) )*dot(vh,n)*ds(1) - dot(n, t(vh, qh, n) )*dot(uh,n)*ds(1) + pen/h*dot(uh,n)*dot(vh, n)*ds(1)

precvarf = inner(Constant(2.)*eta*strain(uh), strain(vh))*dx + Constant(.5)/eta*ph*qh*dx
#weak_bc_prec = -dot(n, tu(uh, n))*dot(vh,n)*ds(1) -dot(n, tu(vh, n) )*dot(uh,n)*ds(1) + pen/h*dot(uh,n)*dot(vh, n)*ds(1)
weak_bc_prec = pen/h*dot(uh,n)*dot(vh, n)*ds(1)

bvarf = dot(rhog,vh) * dx

A = assemble(stokesvarf + weak_bc + robinvarf)
b = assemble(bvarf)
P = assemble(precvarf + weak_bc_prec + robinvarf)

#solver = LUSolver(A)
solver = PETScKrylovSolver("gmres", "amg")
solver.set_operators(A, P)
x = Vector()
A.init_vector(x,0)
n_iter = solver.solve(x,b)
print "Converged in ", n_iter, " iterations."

sol = Function(Xh, x, name = "x")
(u,p) = split(sol)
File("pressure.pvd") << sol.sub(1)
File("velocity.pvd") << sol.sub(0)
plot(u)
plot(p)
interactive()

