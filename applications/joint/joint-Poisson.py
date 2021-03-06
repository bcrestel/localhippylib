import numpy as np
import matplotlib.pyplot as plt
import dolfin as dl
from dolfin import Expression

from hippylib import *
from model_continuous_obs import Poisson
from fenicstools.prior import LaplacianPrior
from fenicstools.regularization import TV, TVPD
from fenicstools.jointregularization import \
SumRegularization, Tikhonovab, VTV, V_TV, V_TVPD
from fenicstools.plotfenics import PlotFenics


if __name__ == "__main__":
    dl.set_log_active(False)
    nx = 64
    ny = 64
    mesh = dl.UnitSquareMesh(nx, ny)
    
    rank = dl.MPI.rank(mesh.mpi_comm())
    nproc = dl.MPI.size(mesh.mpi_comm())
    
    if nproc > 1:
        Random.split(rank, nproc, 1000000, 1)
    
    Vh2 = dl.FunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
    Vh = [Vh2, Vh1, Vh2]
    
    a1true = Expression('log(10' + \
    '- 8*(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4)' + \
    '+ 8*(pow(pow(x[0]-0.25,2)+pow(x[1]-0.5,2),0.5)<0.1)' + \
    '+ 8*((x[0]<=0.8)*(x[0]>=0.7)*(x[1]>=0.45)*(x[1]<=0.55)) )')
    a2true = Expression('log(2' + \
    '+ 8*(pow(pow(x[0]-0.25,2)+pow(x[1]-0.5,2),0.5)<0.1)' + \
    '+ 8*((x[0]<=0.8)*(x[0]>=0.7)*(x[1]>=0.45)*(x[1]<=0.55)) )')
    model1 = Poisson(mesh, Vh, ZeroPrior(Vh[PARAMETER]), atrue=a1true, noiselevel=0.01)
    model2 = Poisson(mesh, Vh, ZeroPrior(Vh[PARAMETER]), atrue=a2true, noiselevel=0.01)
    PltFen = PlotFenics()
    PltFen.set_varname('jointa1')
    PltFen.plot_vtk(model1.at)
    PltFen.set_varname('jointa2')
    PltFen.plot_vtk(model2.at)

    #reg1 = LaplacianPrior({'Vm':Vh[PARAMETER], 'gamma':1e-8, 'beta':1e-8})
    #reg2 = LaplacianPrior({'Vm':Vh[PARAMETER], 'gamma':1e-8, 'beta':1e-8})

    #reg1 = TV({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':1e-8})
    #reg2 = TV({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':1e-8})

    reg1 = TVPD({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':5e-10, 'rescaledradiusdual':1.0})
    reg2 = TVPD({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':5e-10, 'rescaledradiusdual':1.0})

    jointregul = SumRegularization(reg1, reg2, mesh.mpi_comm(), coeff_cg=0.0, coeff_vtv=1.0, \
    parameters_vtv={'eps':1e-3, 'k':5e-9, 'rescaledradiusdual':1.0})
    #jointregul = Tikhonovab({'Vm':Vh[PARAMETER], 'gamma':1e-8, 'beta':1e-8})
    #jointregul = VTV(Vh[PARAMETER], {'k':1e-8, 'eps':1e+1})
    #jointregul = V_TV(Vh[PARAMETER], {'k':1e-8, 'eps':1e-3})
    #jointregul = V_TVPD(Vh[PARAMETER], {'k':5e-9, 'eps':1e-3, 'rescaledradiusdual':1.0})

    jointmodel = JointModel(model1, model2, jointregul)

    solver = ReducedSpaceNewtonCG(jointmodel)
    solver.parameters["rel_tolerance"] = 1e-12
    solver.parameters["abs_tolerance"] = 1e-14
#    solver.parameters["rel_tolerance"] = 1e-10
#    solver.parameters["abs_tolerance"] = 1e-12
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["gda_tolerance"] = 1e-24
    solver.parameters["c_armijo"] = 5e-5
    solver.parameters["max_backtracking_iter"] = 12
    solver.parameters["GN_iter"] = 0
    solver.parameters["max_iter"] = 2000
    solver.parameters["print_level"] = 0
    if rank != 0:
        solver.parameters["print_level"] = -1
    
    InexactCG = 0
    GN = True
    a0 = dl.interpolate(dl.Expression(("0.0","0.0")),jointmodel.Vh[PARAMETER])
    x = solver.solve(a0.vector(), InexactCG, GN)

    x1, x2 = jointmodel.splitvector(x)
    minaf1 = dl.MPI.min(mesh.mpi_comm(), np.amin(x1[PARAMETER].array()))
    maxaf1 = dl.MPI.max(mesh.mpi_comm(), np.amax(x1[PARAMETER].array()))
    md1mis, md1misperc = model1.mediummisfit(x1[PARAMETER])
    minaf2 = dl.MPI.min(mesh.mpi_comm(), np.amin(x2[PARAMETER].array()))
    maxaf2 = dl.MPI.max(mesh.mpi_comm(), np.amax(x2[PARAMETER].array()))
    md2mis, md2misperc = model2.mediummisfit(x2[PARAMETER])
    if rank == 0:
        print 'min(af1)={}, max(af1)={}, med1misft={:e} ({:.1f}%)'.format(\
        minaf1, maxaf1, md1mis, md1misperc)
        print 'min(af2)={}, max(af2)={}, med2misft={:e} ({:.1f}%)'.format(\
        minaf2, maxaf2, md2mis, md2misperc)
        if solver.converged:
            print "\nConverged in ", solver.it, " iterations."
        else:
            print "\nNot Converged"

        print "Termination reason: ", solver.termination_reasons[solver.reason]
        print "Final gradient norm: ", solver.final_grad_norm
        print "Final cost: ", solver.final_cost
    
    PltFen.set_varname('jointsolution1-TVPD+VTVPD')
    PltFen.plot_vtk(vector2Function(x1[PARAMETER], Vh[PARAMETER]))
    PltFen.set_varname('jointsolution2-TVPD+VTVPD')
    PltFen.plot_vtk(vector2Function(x2[PARAMETER], Vh[PARAMETER]))

    if False and nproc == 1:
        xx1 = [vector2Function(x1[i], Vh[i]) for i in range(len(Vh))]
        xx2 = [vector2Function(x2[i], Vh[i]) for i in range(len(Vh))]
        dl.plot(xx1[STATE], title = "State1")
        dl.plot(dl.exp(xx1[PARAMETER]), title = "exp(Parameter1)")
        dl.plot(xx1[ADJOINT], title = "Adjoint1")
        dl.plot(vector2Function(model1.u_o, Vh[STATE]), title = "Observation1")

        dl.plot(xx2[STATE], title = "State2")
        dl.plot(dl.exp(xx2[PARAMETER]), title = "exp(Parameter2)")
        dl.plot(xx2[ADJOINT], title = "Adjoint2")
        dl.plot(vector2Function(model2.u_o, Vh[STATE]), title = "Observation2")
        dl.interactive()


    """
    # Test split and assign
    x = jointmodel.generate_vector("ALL")
    plt.plot(x[STATE].array())
    plt.plot(x[ADJOINT].array())
    plt.plot(x[PARAMETER].array())
    plt.show()

    x1, x2 = jointmodel.splitvector(x, "ALL")
    x1[STATE][:] = 1.
    x1[ADJOINT][:] = 3.
    x1[PARAMETER][:] = 5.
    x2[STATE][:] = 2.
    x2[ADJOINT][:] = 4.
    x2[PARAMETER][:] = 6.
    out = jointmodel.assignvector(x1, x2)
    plt.plot(out[STATE].array(), '-o')
    plt.plot(out[ADJOINT].array(), '-o')
    plt.plot(out[PARAMETER].array(), '-o')
    plt.show()
    """
