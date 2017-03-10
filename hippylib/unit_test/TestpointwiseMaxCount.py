import dolfin as dl
from hippylib import *

mesh = dl.UnitSquareMesh(3,3)
rank = dl.MPI.rank(mesh.mpi_comm())
print 'rank={}'.format(rank)

V = dl.FunctionSpace(mesh,'CG',1)
test = dl.TestFunction(V)
trial = dl.TrialFunction(V)
M = dl.assemble(dl.inner(test,trial)*dl.dx)
u = dl.interpolate(dl.Expression("1+sin(2*x[0]*pi)*sin(2*x[1]*pi)"), V).vector()
v = dl.Vector()
M.init_vector(v,0)
print 'rank={}'.format(rank)

#print u.array()
count = pointwiseMaxCount(v, u, 0.3)
#print u.array()
print 'rank={}, count={}'.format(rank, count)
print v.array()
