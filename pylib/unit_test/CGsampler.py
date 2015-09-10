import dolfin as dl
import sys
sys.path.append( "../../" )
from pylib import *
import numpy as np

nx = 32
ny = 32
mesh = dl.UnitSquareMesh(nx, ny)
Vh = dl.FunctionSpace(mesh, 'Lagrange', 1)
Prior = LaplacianPrior(Vh, 1.,100.)

nsamples = 1000

s = dl.Function(Vh, name = "sample")
noise = dl.Vector()
Prior.init_vector(noise,"noise")
size = len(noise.array())

fid = dl.File("results_cg/samples.pvd")

for i in range(0, nsamples ):
    noise.set_local( np.random.randn( size ) )
    Prior.sample(noise, s.vector())
    fid << s


