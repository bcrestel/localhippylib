import numpy as np
import matplotlib.pyplot as plt
import dolfin as dl
from dolfin import Expression

from hippylib import *
from model_pointwise_obs import Poisson

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
    
    a1true = dl.Expression('log(10 - ' + \
    '(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4) * (' + \
    '4*(x[0]<=0.5) + 8*(x[0]>0.5) ))')
    a2true = dl.Expression('log(10 - ' + \
    '(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4) * (' + \
    '8*(x[0]<=0.5) + 4*(x[0]>0.5) ))')

    nbobsperdir=50
    targets1 = np.array([ [float(i)/(nbobsperdir+1), float(j)/(nbobsperdir+1)] \
    for i in range((nbobsperdir+2)/2, nbobsperdir+1) \
    for j in range((nbobsperdir+2)/2, nbobsperdir+1)])
    targets2 = np.array([ [float(i)/(nbobsperdir+1), float(j)/(nbobsperdir+1)] \
    for i in range(1, nbobsperdir+1) for j in range(1, nbobsperdir+1)])

    model1 = Poisson(mesh, Vh, a1true, targets1, ZeroPrior(Vh[PARAMETER]), 0.02)
    model2 = Poisson(mesh, Vh, a2true, targets2, ZeroPrior(Vh[PARAMETER]), 0.02)
    PltFen = PlotFenics()
    PltFen.set_varname('jointa1')
    PltFen.plot_vtk(model1.at)
    PltFen.set_varname('jointa2')
    PltFen.plot_vtk(model2.at)

    #reg1 = LaplacianPrior({'Vm':Vh[PARAMETER], 'gamma':1e-8, 'beta':1e-8})
    #reg2 = LaplacianPrior({'Vm':Vh[PARAMETER], 'gamma':1e-8, 'beta':1e-8})

    #reg1 = TV({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':1e-8})
    #reg2 = TV({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':1e-8})

    reg1 = TVPD({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':3e-7, 'rescaledradiusdual':1.0})
    reg2 = TVPD({'Vm':Vh[PARAMETER], 'eps':1e-3, 'k':4e-7, 'rescaledradiusdual':1.0})

    jointregul = SumRegularization(reg1, reg2, mesh.mpi_comm(), coeff_cg=1e-4, coeff_vtv=0.0, \
    parameters_vtv={'eps':1e-3, 'k':5e-9, 'rescaledradiusdual':1.0})
    #jointregul = Tikhonovab({'Vm':Vh[PARAMETER], 'gamma':1e-8, 'beta':1e-8})
    #jointregul = VTV(Vh[PARAMETER], {'k':1e-8, 'eps':1e+1})
    #jointregul = V_TV(Vh[PARAMETER], {'k':4e-7, 'eps':1e-3})
    #jointregul = V_TVPD(Vh[PARAMETER], {'k':4e-7, 'eps':1e-7, 'rescaledradiusdual':1.0})

    ##### Modify this #####
    plot_suffix = 'TVPD+1e-4CG'
    #######################

    jointmodel = JointModel(model1, model2, jointregul)

    solver = ReducedSpaceNewtonCG(jointmodel)
    solver.parameters["rel_tolerance"] = 1e-12
    solver.parameters["abs_tolerance"] = 1e-14
#    solver.parameters["rel_tolerance"] = 1e-10
#    solver.parameters["abs_tolerance"] = 1e-12
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["gda_tolerance"] = 1e-24
    solver.parameters["c_armijo"] = 5e-5
    solver.parameters["max_backtracking_iter"] = 20 # !!! very large
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
    
    PltFen.set_varname('jointsolution1-' + plot_suffix)
    PltFen.plot_vtk(vector2Function(x1[PARAMETER], Vh[PARAMETER]))
    PltFen.set_varname('jointsolution2-' + plot_suffix)
    PltFen.plot_vtk(vector2Function(x2[PARAMETER], Vh[PARAMETER]))
