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
import sys
import dolfin as dl
import math
import numpy as np
import matplotlib.pyplot as plt

import sys, os
from hippylib import Random, ReducedSpaceNewtonCG, amg_method, Model,\
STATE, ADJOINT, PARAMETER, vector2Function

from fenicstools.regularization import TVPD
from fenicstools.plotfenics import PlotFenics

from definePDE import pdes
from definemisfit import defmisfit
from targetmedium import targetmediumparameters, initmediumparameters

dl.set_log_active(False)


PLOT = True
            
sep = "\n"+"#"*80+"\n"
ndim = 2
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

ndofs = [Vh[STATE].dim(), Vh[PARAMETER].dim(), Vh[ADJOINT].dim()]
if rank == 0:
    print sep, "Set up the mesh and finite element spaces", sep
    print "Number of dofs: STATE={0}, PARAMETER={1}, ADJOINT={2}".format(*ndofs)

# Target medium parameters
atrue,_,_,_,_ = targetmediumparameters(Vh[PARAMETER], 1.0)
    
# Define PDE
pde = pdes(Vh, STATE, amg_method)

# Define misfit functions
misfit, targets = defmisfit(Vh, STATE)

# Generate synthetic observations
rel_noise_level = 0.02
utrue = pde.generate_state()
x = [utrue, atrue.vector(), None]
minatrue = dl.MPI.min(mesh.mpi_comm(), np.amin(atrue.vector().array()))
maxatrue = dl.MPI.max(mesh.mpi_comm(), np.amax(atrue.vector().array()))
if rank == 0:
    print 'min(atrue)={}, max(atrue)={}'.format(minatrue, maxatrue)
pde.solveFwd(x[STATE], x, 1e-9)
noise_level = rel_noise_level * x[STATE].norm("l2") / np.sqrt(Vh[PARAMETER].dim())
Random.normal(x[STATE], noise_level, False)
misfit.B.mult(x[STATE], misfit.d)
#misfit.noise_variance = np.sqrt(targets.shape[0])   # hack to compare both models
misfit.noise_variance = 1.0

# Regularization
prior = TVPD({'Vm':Vh[PARAMETER], 'k':1e-3, 'eps':1e-5, 'print':(not rank)})

if PLOT:
    PltFen = PlotFenics(mesh.mpi_comm(), os.path.splitext(sys.argv[0])[0] + '/Plots')
    suffix = '-c1-k' + str(prior.parameters['k']) + \
    '-e' + str(prior.parameters['eps'])
    PltFen.set_varname('atrue')
    PltFen.plot_vtk(atrue)

model = Model(pde, prior, misfit, atrue.vector())
model.solveFwd(x[STATE], x, 1e-9)
c, r, m = model.cost(x)
if rank == 0:
    print 'Cost @ MAP: cost={}, misfit={}, reg={}'.format(c, m, r)

if rank == 0:
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
if rank != 0:
    solver.parameters["print_level"] = -1

#a0,_,_,_,_ = initmediumparameters(Vh[PARAMETER], 1.0)
a0 = dl.interpolate(dl.Constant('0.0'), Vh[PARAMETER])
x = solver.solve(a0.vector(), InexactCG=1, GN=False, bounds_xPARAM=[-10, 15])
#x = solver.solve(a0.vector(), InexactCG=1, GN=False, bounds_xPARAM=[1e-4, 1.0])

minaf = dl.MPI.min(mesh.mpi_comm(), np.amin(x[PARAMETER].array()))
maxaf = dl.MPI.max(mesh.mpi_comm(), np.amax(x[PARAMETER].array()))
mdmis, mdmisperc = model.mediummisfit(x[PARAMETER])
if rank == 0:
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
    PltFen.set_varname('aMAP')
    PltFen.plot_vtk(vector2Function(x[PARAMETER], Vh[PARAMETER]))
