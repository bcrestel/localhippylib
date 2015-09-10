import dolfin as dl
import sys
sys.path.append( "../" )
from pylib import *
import numpy as np

dl.set_log_active(False)
nx = 4
ny = 4
nz = 1
mesh = dl.UnitCubeMesh(nx, ny, nz)

for i in range(2):
    cell_markers = dl.CellFunction("bool", mesh)
    cell_markers.set_all(False)
    for cell in dl.cells(mesh):
        if cell.midpoint()[1] < .75 and cell.midpoint()[1] > .25 and cell.midpoint()[0] > .2 and cell.midpoint()[0] < .5:
            cell_markers[cell] = True
            
    mesh = dl.refine(mesh, cell_markers)
    
order = 2

Vh = dl.FunctionSpace(mesh, 'Lagrange', order)
X1h = dl.FunctionSpace(mesh, 'Quadrature', 2*order)
Xh = dl.MixedFunctionSpace( [X1h, X1h, X1h, X1h] )
uh = dl.TrialFunction(Vh)
vh = dl.TestFunction(Vh)
xh = dl.TrialFunction(Xh)
zh = dl.TestFunction(Xh)
z1h, z2h, z3h, z4h = dl.split(zh)

dl.plot(mesh, interactive=True)

A = dl.assemble(uh*vh*dl.dx + dl.inner(dl.nabla_grad(uh), dl.nabla_grad(vh))*dl.dx )
W = dl.assemble( dl.inner(xh,zh)*dl.dx )
L = dl.assemble( uh.dx(0)*z1h*dl.dx + uh.dx(1)*z2h*dl.dx + uh.dx(2)*z3h*dl.dx + uh*z4h*dl.dx )

ones = dl.Vector()
W.init_vector(ones,0)
ones.set_local(np.ones(ones.size()))
W1 = W*ones
invW1 = dl.Vector()
W.init_vector(invW1,0)
invW1.set_local(1./W1.array())

Winv = dl.assemble( dl.inner(xh,zh)*dl.dx )
Winv.zero()
Winv.set_diagonal(invW1)

randx = dl.Vector()
A.init_vector(randx,0)
randx.set_local(np.random.rand(randx.size() ) )

for x in [dl.interpolate(dl.Constant(1), Vh).vector(),
          dl.interpolate(dl.Expression("x[0]"), Vh).vector(),
          dl.interpolate(dl.Expression("x[0]*x[1]"), Vh).vector(),
          dl.interpolate(dl.Expression("x[0]*x[0]"), Vh).vector(),
          dl.interpolate(dl.Expression("x[0]*x[0]*x[1]"), Vh).vector(),
          randx]:
    y = A*x
    Lx = L*x
    help = Winv * Lx
    y2 = dl.Vector()
    L.transpmult(help, y2)
    print (y - y2).norm("linf")/y.norm("linf")
