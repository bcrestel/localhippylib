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
import os
    
abspath = os.path.dirname( os.path.abspath(__file__) )
sdir = os.path.join(abspath,"AssemblePointwiseObservation")
header_file = open(os.path.join(sdir,"AssemblePointwiseObservation.h"), "r")
code = header_file.read()
header_file.close()
cpp_sources = ["AssemblePointwiseObservation.cpp"]
cpp_module = dl.compile_extension_module(
code=code, source_directory=sdir, sources=cpp_sources,
include_dirs=[".",  sdir])

def assemblePointwiseObservation(Vh, targets):
    """
    Assemble the pointwise observation matrix:
    Input
    - Vh: FEniCS finite element space
    - targets: observation points (numpy array)
     
    Note: This Function will not work in parallel!!!
    """
    #Ensure that PetscInitialize is called
    dummy = dl.assemble( dl.inner(dl.TrialFunction(Vh), dl.TestFunction(Vh))*dl.dx )
    #Call the cpp module to compute the pointwise observation matrix
    tmp = cpp_module.PointwiseObservation(Vh,targets.flatten())
    #return the matrix
    return tmp.GetMatrix()

def exportPointwiseObservation(Vh, B, data, fname, varname="observation"):
    """
    This function write a VTK PolyData file to visualize pointwise data.
    Inputs:
    - B: observation operator
    - mesh: mesh
    - data: dl.Vector containing the data
    - fname:   filename for the vtk file to export
    - varname: name of the variable for the vtk file
    """
    mesh = Vh.mesh()
    if mesh.geometry().dim() == 1:
        xyz_fun = [dl.Expression("x[0]", degree=1)]
    elif mesh.geometry().dim() == 1:
        xyz_fun = [dl.Expression("x[0]", degree=1), dl.Expression("x[1]", degree=1)]
    else:
        xyz_fun = [dl.Expression("x[0]", degree=1), dl.Expression("x[1]", degree=1), dl.Expression("x[2]", degree=1)]

    xyz = [B*dl.interpolate(fun, Vh).vector() for fun in xyz_fun]
    
    data_on_pzero = data.gather_on_zero()
    xyz_on_pzero = np.zeros((data_on_pzero.shape[0], 3))
    for i in range(len(xyz)):
        xyz_on_pzero[:,i] = xyz[i].gather_on_zero()
    
    #check if I'm rank 0, walkaroud to avoid mpi4py
    if B.local_range(0)[0] == 0:
        write_vtk(xyz_on_pzero, data_on_pzero, fname, varname)
    

def write_vtk(points, data, fname, varname="observation"):
    """
    This function write a VTK PolyData file to visualize pointwise data.
    Inputs:
    - points:  locations of the points 
               (numpy array of size number of points by space dimension)
    - data:    pointwise values
               (numpy array of size number of points)
    - fname:   filename for the vtk file to export
    - varname: name of the variable for the vtk file
    """
    ndim = points.shape[1]
    npoints = points.shape[0]
    
    assert npoints == data.shape[0]
    assert ndim == 3
    
    f = open(fname, 'w')
    f.write('<VTKFile version="0.1" byte_order="LittleEndian" type="PolyData">\n')
    f.write('<PolyData>\n')
    f.write('<Piece NumberOfPoints="{0:d}" NumberOfVerts="{1:d}">\n'.format(npoints,npoints))
    f.write('<Points>\n')
    f.write('<DataArray NumberOfComponents="{0:d}" format="ascii" type="Float32">\n'.format(ndim))
    f.write( np.array_str(points, precision=9, suppress_small=True).replace("[", "").replace("]", "") )
    f.write('\n</DataArray>\n')
    f.write('</Points>\n')
    f.write('<Verts>\n')
    f.write('\n<DataArray type="Int32" Name="connectivity" format="ascii">\n')
    f.write(np.array_str( np.arange(0,npoints) ).replace("[", "").replace("]", "") )
    f.write('</DataArray>\n')
    f.write('<DataArray type="Int32" Name="offsets" format="ascii">\n')
    f.write(np.array_str( np.arange(1,npoints+1) ).replace("[", "").replace("]", "") )
    f.write('</DataArray>\n')
    f.write('</Verts>\n')
    f.write('<PointData Scalars="{}">\n'.format(varname))
    f.write('<DataArray format="ascii" type="Float32" Name="{}">\n'.format(varname))
    f.write(np.array_str( data, precision=9, suppress_small=True ).replace("[", "").replace("]", "") )
    f.write('\n</DataArray>\n</PointData>\n<CellData>\n</CellData>\n</Piece>\n</PolyData>\n</VTKFile>')
    f.close()