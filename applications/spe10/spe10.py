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


from dolfin.cpp.mesh import RectangleMesh

"""
This script reads the SPE10 benchmark permeability field from file and it
creates a dolfin::Expression to evaluate the permeability field at a desired
point in the the domain.
"""

from dolfin import *
import math

class SPE10Coefficient3D(Expression):
    def setFilename(self, filename):
        # number of cells in each dimension
        self.N = [60, 220, 85]
        # lengths of the domain in ft
        self.L = [1200, 2200, 170]
        # h
        self.h = [float(l)/float(n) for (l,n) in zip(self.L, self.N)]
        
        # This opens a handle to your file, in 'r' read mode
        file_handle = open(filename, 'r')
        # Read in all the lines of your file into a list of lines
        lines_list = file_handle.readlines()
        # extract all the permeability values from each line
        self.data = [ float(val) for line in lines_list for val in line.split()]
        
        
    def eval(self, value, x):
        #find to which cell i,i,k does x belongs
        i = int( math.floor( x[0]/self.h[0] ) )
        j = int( math.floor( x[1]/self.h[1] ) )
        k = self.N[2] - 1 - int( math.floor( x[2]/self.h[2] ) )
        
        # ensure that the point is valid (i.e. we are not extrapolating)
        assert( i < self.N[0] )
        assert( i >= 0 )
        assert( j < self.N[1] )
        assert( j >= 0 )
        assert( k < self.N[2] )
        assert( k >= 0 )
        
        #compute the index in the self.data list (x is the fastest running coordinate, then y, than z)
        index = self.N[0]*self.N[1]*k + self.N[0]*j + i   
        value[0] = self.data[index]
        value[1] = self.data[index + self.N[0]*self.N[1]*self.N[2] ]
        value[2] = self.data[index + 2*self.N[0]*self.N[1]*self.N[2] ]

    def value_shape(self):
        return (3,)
    
class SPE10Coefficient2D(Expression):
    def setFilename(self, filename, zslice=0.):
        # number of cells in each dimension
        self.N = [60, 220, 85]
        # lengths of the domain in ft
        self.L = [1200, 2200, 170]
        # h
        self.h = [float(l)/float(n) for (l,n) in zip(self.L, self.N)]
        
        # This opens a handle to your file, in 'r' read mode
        file_handle = open(filename, 'r')
        # Read in all the lines of your file into a list of lines
        lines_list = file_handle.readlines()
        # extract all the permeability values from each line
        self.data = [ float(val) for line in lines_list for val in line.split()]
        self.zslice = zslice
        
    def eval(self, value, x):
        #find to which cell i,i,k does x belongs
        i = int( math.floor( x[0]/self.h[0] ) )
        j = int( math.floor( x[1]/self.h[1] ) )
        k = self.N[2] - 1 - int( math.floor( self.zslice/self.h[2] ) )
        
        # ensure that the point is valid (i.e. we are not extrapolating)
        assert( i < self.N[0] )
        assert( i >= 0 )
        assert( j < self.N[1] )
        assert( j >= 0 )
        assert( k < self.N[2] )
        assert( k >= 0 )
        
        #compute the index in the self.data list (x is the fastest running coordinate, then y, than z)
        index = self.N[0]*self.N[1]*k + self.N[0]*j + i   
        value[0] = self.data[index]


    
if __name__ == "__main__":
    # define the mesh
    twod = True
    
    if twod:
        mesh = RectangleMesh(0., 0., 1200, 2200, 60, 220)
        
        perm = SPE10Coefficient2D()
        perm.setFilename("por_perm_case2a/spe_perm.dat")
        
        #define a finite element space to interpolate the coefficient
        Vh = FunctionSpace(mesh, 'DG', 0)
        perm_h = interpolate(perm, Vh)
        #    plot(ch, interactive=True)
        plot(ln(perm_h), interactive=True)
        
        
    else:
        mesh = BoxMesh(0., 0., 0., 1200, 2200, 170, 60, 220, 85)
        
        perm = SPE10Coefficient3D()
        perm.setFilename("por_perm_case2a/spe_perm.dat") 
        
        #define a finite element space to interpolate the coefficient
        Vh = VectorFunctionSpace(mesh, 'DG', 0)
        perm_h = interpolate(perm, Vh)
        #    plot(ch, interactive=True)
        plot(ln(perm_h.sub(0)), interactive=True)
        plot(ln(perm_h.sub(1)), interactive=True)
        plot(ln(perm_h.sub(2)), interactive=True)
