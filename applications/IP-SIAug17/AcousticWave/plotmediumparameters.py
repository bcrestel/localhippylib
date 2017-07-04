"""
Plot selected medium parameters
"""
import sys
from os.path import splitext, isdir
from shutil import rmtree
import dolfin as dl
from dolfin import MPI
from fenicstools.plotfenics import PlotFenics
from targetmedium import \
targetmediumparameters, initmediumparameters, loadparameters

LARGE = False
Nxy, Dt, fpeak,_,_,_,tf = loadparameters(LARGE)

mesh = dl.UnitSquareMesh(Nxy, Nxy)
Vl = dl.FunctionSpace(mesh, 'Lagrange', 1)

mpicomm = mesh.mpi_comm()
mpirank = MPI.rank(mpicomm)
filename, ext = splitext(sys.argv[0])
if mpirank == 0: 
    if isdir(filename + '/'):   rmtree(filename + '/')
MPI.barrier(mpicomm)
myplot = PlotFenics(mpicomm, filename)

at, bt = targetmediumparameters(Vl, 1.0, myplot)
a0, b0 = initmediumparameters(Vl, 1.0, myplot)
at, bt = targetmediumparameters(Vl, 0.0, myplot)
a0, b0 = initmediumparameters(Vl, 0.0, myplot)
