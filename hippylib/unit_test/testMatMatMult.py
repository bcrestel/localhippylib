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
import sys
sys.path.append( "../../" )
from hippylib import *
import numpy as np

dl.set_log_active(False)
nx = 8
ny = 8
mesh = dl.UnitSquareMesh(nx, ny)

Vh = dl.FunctionSpace(mesh, 'Lagrange', 1)

data = np.ones(Vh.dim(), dtype=np.float64)

uh = dl.TrialFunction(Vh)
vh = dl.TestFunction(Vh)

A = dl.assemble(uh*vh*dl.dx)
B = dl.assemble(dl.inner(dl.nabla_grad(uh), dl.nabla_grad(vh))*dl.dx)

C = MatMatMult(A,B)
D = MatPtAP(A,B)


x = dl.Vector()
y = dl.Vector()
xA = dl.Vector()
yA = dl.Vector()
xD = dl.Vector()
yD = dl.Vector()

A.init_vector(yA,0)
A.init_vector(xA,1)

C.init_vector(y,0)
C.init_vector(x,1)

D.init_vector(yD,0)
D.init_vector(xD,1)

print "y in range of A: y.set_local(data)..."
yA.set_local(data)
print "... success"

print "x in domain of A: x.set_local(data)..."
xA.set_local(data)
print "... success"

print "x in domain of C: x.set_local(data)..."
x.set_local(data)
print "... success"

print "x in domain of C: x.set_local(data)..."
x.set_local(data)
print "... success"

print "y in range of C: y.set_local(data)..."
y.set_local(data)
print "... success"

print "x in domain of C: x.set_local(data)..."
x.set_local(data)
print "... success"

print "y in range of D: y.set_local(data)..."
yD.set_local(data)
print "... success"

print "x in domain of D: x.set_local(data)..."
xD.set_local(data)
print "... success"
