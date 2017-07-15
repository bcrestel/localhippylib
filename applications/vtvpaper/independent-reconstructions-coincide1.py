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
import math
import numpy as np
import matplotlib.pyplot as plt

import sys, os
from hippylib import *

from fenicstools.regularization import TVPD
from fenicstools.plotfenics import PlotFenics

from targetmedium_coincide1 import targetmedium
from definePDE import pdes
from definemisfit import defmisfit



PLOT = True
            
if __name__ == "__main__":
    dl.set_log_active(False)
    sep = "\n"+"#"*80+"\n"
    ndim = 2
    nx = 64
    ny = 64

    ######################################
    try:
        SELECTMODEL = int(sys.argv[1])
        k = float(sys.argv[2])
        eps = float(sys.argv[3])
    except:
        SELECTMODEL = 1
        k = 3e-7
        eps = 1e-3
    ######################################
    
    mesh = dl.UnitSquareMesh(nx, ny)

    mpicomm = mesh.mpi_comm()
    mpirank = dl.MPI.rank(mpicomm)
    mpisize = dl.MPI.size(mpicomm)

    PRINT = (mpirank == 0)
    
    if mpisize > 1:
        Random.split(mpirank, mpisize, 1000000, 1)
        
    Vh2 = dl.FunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
    Vh = [Vh2, Vh1, Vh2]
    
    ndofs = [Vh[STATE].dim(), Vh[PARAMETER].dim(), Vh[ADJOINT].dim()]
    if PRINT:
        print sep, "Set up the mesh and finite element spaces", sep
        print "Number of dofs: STATE={0}, PARAMETER={1}, ADJOINT={2}".format(*ndofs)
    
    # Target medium parameters
    a1true, a2true = targetmedium(Vh, PARAMETER)
        
    # Define PDE
    pde1, pde2 = pdes(Vh, STATE, amg_method)
 
    # Define misfit functions
    misfit1, misfit2, targets1, targets2 = defmisfit(Vh, STATE)

    # Generate synthetic observations
    rel_noise_level = 0.02
    utrue1 = pde1.generate_state()
    utrue2 = pde2.generate_state()
    for misfit, atrue, targets, utrue, pde in zip([misfit1, misfit2], \
    [a1true, a2true], [targets1, targets2], [utrue1, utrue2], [pde1, pde2]):
        x = [utrue, atrue.vector(), None]
        minatrue = dl.MPI.min(mesh.mpi_comm(), np.amin(atrue.vector().array()))
        maxatrue = dl.MPI.max(mesh.mpi_comm(), np.amax(atrue.vector().array()))
        if PRINT:
            print 'min(atrue)={}, max(atrue)={}'.format(minatrue, maxatrue)
        pde.solveFwd(x[STATE], x, 1e-9)
        noise_level = rel_noise_level * x[STATE].norm("l2") / np.sqrt(Vh[PARAMETER].dim())
        Random.normal(x[STATE], noise_level, False)
        misfit.B.mult(x[STATE], misfit.d)
        misfit.noise_variance = np.sqrt(targets.shape[0])   # hack to compare both models
    
    # Regularization
    prior = TVPD({'Vm':Vh[PARAMETER], 'k':k, 'eps':eps, 'print':PRINT})

    if PLOT:
        PltFen = PlotFenics(mesh.mpi_comm(), os.path.splitext(sys.argv[0])[0] + '/plots')
        filename = 'model' + str(SELECTMODEL) + '-TVPD_k' + str(k) + '_eps' + str(eps)
        PltFen.set_varname('a1')
        PltFen.plot_vtk(a1true)
        PltFen.set_varname('a2')
        PltFen.plot_vtk(a2true)

    if SELECTMODEL == 1:
        atrue = a1true
        pde = pde1
        misfit = misfit1
    else:
        atrue = a2true
        pde = pde2
        misfit = misfit2

    model = Model(pde, prior, misfit, atrue.vector())
    x[STATE].zero()
    c, r, m = model.cost(x)
    if PRINT:
        print 'Cost @ MAP: cost={}, misfit={}, reg={}'.format(c, m, r)
    
    if PRINT:
        print sep, "Find the MAP point", sep
    solver = ReducedSpaceNewtonCG(model)
    solver.parameters["rel_tolerance"] = 1e-12
    solver.parameters["abs_tolerance"] = 1e-14
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["gda_tolerance"] = 1e-24
    solver.parameters["c_armijo"] = 5e-5
    solver.parameters["max_backtracking_iter"] = 20
    solver.parameters["GN_iter"] = 5
    solver.parameters["max_iter"] = 2000
    solver.parameters["print_level"] = 0
    if not PRINT:
        solver.parameters["print_level"] = -1
    
    a0 = dl.interpolate(dl.Constant('0.0'), Vh[PARAMETER])
    x = solver.solve(a0.vector(), InexactCG=1, GN=False, bounds_xPARAM=[-10., 25.])
    
    minaf = x[PARAMETER].min()
    maxaf = x[PARAMETER].max()
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
        PltFen.set_varname(filename)
        PltFen.plot_vtk(vector2Function(x[PARAMETER], Vh[PARAMETER]))
