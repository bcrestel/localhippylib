import matplotlib.pyplot as plt
import dolfin as dl

from hippylib import *
from model_continuous_obs import Poisson
from fenicstools.prior import LaplacianPrior
from fenicstools.jointregularization import Tikhonovab


if __name__ == "__main__":
    dl.set_log_active(False)
    nx = 4
    ny = 4
    mesh = dl.UnitSquareMesh(nx, ny)
    
    rank = dl.MPI.rank(mesh.mpi_comm())
    nproc = dl.MPI.size(mesh.mpi_comm())
    
    if nproc > 1:
        Random.split(rank, nproc, 1000000, 1)
    
    Vh2 = dl.FunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
    Vh = [Vh2, Vh1, Vh2]
    
    model1 = Poisson(mesh, Vh, ZeroPrior(Vh[PARAMETER]))
    model2 = Poisson(mesh, Vh, ZeroPrior(Vh[PARAMETER]))

    jointregul = Tikhonovab({'Vm':Vh[PARAMETER], 'gamma':1e-7. 'beta':1e-8})

    jointmodel = JointModel(model1, model2, jointregul)

    a0 = interpolate(Expression("0.0"),Vh[PARAMETER])
    solver = ReducedSpaceNewtonCG(model)
    solver.parameters["rel_tolerance"] = 1e-10
    solver.parameters["abs_tolerance"] = 1e-12
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["c_armijo"] = 5e-5
    solver.parameters["max_backtracking_iter"] = 12
    solver.parameters["GN_iter"] = 5
    solver.parameters["max_iter"] = 2000
    if rank != 0:
        solver.parameters["print_level"] = -1
    
    InexactCG = 1
    GN = True
    #x = solver.solve(a0.vector(), InexactCG, GN)



#    # Test split and assign
#    x = jointmodel.generate_vector("ALL")
#    plt.plot(x[STATE].array())
#    plt.plot(x[ADJOINT].array())
#    plt.plot(x[PARAMETER].array())
#    plt.show()
#
#    x1, x2 = jointmodel.splitvector(x, "ALL")
#    x1[STATE][:] = 1.
#    x1[ADJOINT][:] = 3.
#    x1[PARAMETER][:] = 5.
#    x2[STATE][:] = 2.
#    x2[ADJOINT][:] = 4.
#    x2[PARAMETER][:] = 6.
#    out = jointmodel.assignvector(x1, x2)
#    plt.plot(out[STATE].array(), '-o')
#    plt.plot(out[ADJOINT].array(), '-o')
#    plt.plot(out[PARAMETER].array(), '-o')
#    plt.show()
