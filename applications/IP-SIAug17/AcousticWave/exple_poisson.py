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
# Software Foundation) version 3.0 dated June 2007.
import dolfin as dl
import numpy as np

import sys, os
from hippylib import ReducedSpaceNewtonCG, amg_method, Model,\
STATE, ADJOINT, PARAMETER, vector2Function

from fenicstools.regularization import TVPD
from fenicstools.plotfenics import PlotFenics

from definePDE import pdes
from definemisfit import defmisfit
from targetmedium import targetmediumparameters, initmediumparameters

dl.set_log_active(False)



def model_poisson(Vh, prior, PRINT=False):
    # Target medium parameters
    atrue,_ = targetmediumparameters(Vh[PARAMETER], 0.0)

    # Define PDE
    pde = pdes(Vh, STATE, amg_method)

    # Define misfit functions
    misfit,_ = defmisfit(Vh, STATE)

    # Generate synthetic observations
    rel_noise_level = 0.01
    utrue = pde.generate_state()
    rnd_v = pde.generate_state()
    x = [utrue, atrue.vector(), None]
    mpicomm = Vh[PARAMETER].mesh().mpi_comm()
    minatrue = dl.MPI.min(mpicomm, np.amin(atrue.vector().array()))
    maxatrue = dl.MPI.max(mpicomm, np.amax(atrue.vector().array()))
    if PRINT:
        print '[poisson] min(atrue)={}, max(atrue)={}'.format(minatrue, maxatrue)
    pde.solveFwd(x[STATE], x, 1e-9)
    noise_level = rel_noise_level * x[STATE].norm("l2") / np.sqrt(Vh[PARAMETER].dim())
    np.random.seed(1111)
    rnd = np.random.randn(x[STATE].local_size())
    rnd_v[:] = rnd
    x[STATE].axpy(noise_level, rnd_v)
    misfit.B.mult(x[STATE], misfit.d)
    misfit.noise_variance = 5e3 # hack to compare elliptic and acoustic
    #misfit.noise_variance = 1.0

    model = Model(pde, prior, misfit, atrue.vector())
    model.solveFwd(x[STATE], x, 1e-9)
    c, r, m = model.cost(x)
    if PRINT:
        print '[poisson] Cost @ MAP: cost={}, misfit={}, reg={}'.format(c, m, r)

    return model




if __name__ == "__main__":

    # Command-line argument
    try:
        k = float(sys.argv[1])
        eps = float(sys.argv[2])
    except:
        k = 2e-8
        eps = 1e-3
    #######################


    PLOT = True
                
    sep = "\n"+"#"*80+"\n"

    ndim = 2

    nx = 20
    ny = 20
    mesh = dl.UnitSquareMesh(nx, ny)

    rank = dl.MPI.rank(mesh.mpi_comm())
    nproc = dl.MPI.size(mesh.mpi_comm())
    PRINT = (rank == 0)

    Vh2 = dl.FunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
    Vh = [Vh2, Vh1, Vh2]

    ndofs = [Vh[STATE].dim(), Vh[PARAMETER].dim(), Vh[ADJOINT].dim()]
    if PRINT:
        print sep, "Set up the mesh and finite element spaces", sep
        print "Number of dofs: STATE={0}, PARAMETER={1}, ADJOINT={2}".format(*ndofs)

    # Regularization
    prior = TVPD({'Vm':Vh[PARAMETER], 'k':k, 'eps':eps, 'print':(not rank)})

    model = model_poisson(Vh, prior, PRINT)

    if PLOT:
        PltFen = PlotFenics(mesh.mpi_comm(), os.path.splitext(sys.argv[0])[0] + '/plots')
        suffix = '-c1-k' + str(prior.parameters['k']) + \
        '-e' + str(prior.parameters['eps'])
        PltFen.set_varname('atrue')
        PltFen.plot_vtk(vector2Function(model.atrue, Vh[PARAMETER]))

    if PRINT:
        print sep, "Find the MAP point", sep
    solver = ReducedSpaceNewtonCG(model)
    solver.parameters["rel_tolerance"] = 1e-12
    solver.parameters["abs_tolerance"] = 1e-14
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["gda_tolerance"] = 1e-24
    solver.parameters["c_armijo"] = 5e-5
    solver.parameters["max_backtracking_iter"] = 20
    solver.parameters["GN_iter"] = 10
    solver.parameters["max_iter"] = 2000
    solver.parameters["print_level"] = 0
    if not PRINT:
        solver.parameters["print_level"] = -1

    a0,_ = initmediumparameters(Vh[PARAMETER], 0.0)
    x = solver.solve(a0.vector(), InexactCG=1, GN=False, bounds_xPARAM=[1e-8, 100.0])

    minaf = dl.MPI.min(mesh.mpi_comm(), np.amin(x[PARAMETER].array()))
    maxaf = dl.MPI.max(mesh.mpi_comm(), np.amax(x[PARAMETER].array()))
    mdmis, mdmisperc = model.mediummisfit(x[PARAMETER])
    if PRINT:
        print 'min(af)={}, max(af)={}, medmisft={:e} ({:.1f}%)'.format(\
        minaf, maxaf, mdmis, mdmisperc)
        if solver.converged:
            print "\nConverged in ", solver.it, " iterations."
        else:
            print "\nNot Converged"

        print "Termination reason: ", solver.termination_reasons[solver.reason]
        print "Final gradient norm: ", solver.final_grad_norm
        print "Final cost: ", solver.final_cost

    # Plot reconstruction
    if PLOT:
        PltFen.set_varname('poisson-MAP_k' + str(k) + '_e' + str(eps))
        PltFen.plot_vtk(vector2Function(x[PARAMETER], Vh[PARAMETER]))
