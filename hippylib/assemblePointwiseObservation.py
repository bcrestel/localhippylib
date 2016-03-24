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
# Software Foundation) version 2.1 dated February 1999.

import dolfin as dl
import os
    
abspath = os.path.dirname( os.path.abspath(__file__) )
sdir = os.path.join(abspath,"AssemblePointwiseObservation")
header_file = open(os.path.join(sdir,"AssemblePointwiseObservation.h"), "r")
code = header_file.read()
header_file.close()
#check the dolfin version to decide which cpp to include
if dl.dolfin_version()[2] == "4":
    cpp_sources = ["AssemblePointwiseObservation_v14.cpp"]
elif dl.dolfin_version()[2] == "5":
    cpp_sources = ["AssemblePointwiseObservation_v15.cpp"]
elif dl.dolfin_version()[2] == "6":
    cpp_sources = ["AssemblePointwiseObservation_v16.cpp"]
else:
    raise Exception("Dolfin Version")

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
    
    

#    from petsc4py import PETSc
#    import numpy as np
#
#    def assemblePointwiseObservation(Vh, targets):
#        """
#        Assemble the pointwise observation matrix:
#        Input
#        - Vh: FEniCS finite element space
#        - targets: observation points (numpy array)
#         
#        Note: This Function will not work in parallel!!!
#        """
#        ntargets, dim = targets.shape
#        mesh = Vh.mesh()
#        coords = mesh.coordinates()
#        cells = mesh.cells()
#        dolfin_element = Vh.dolfin_element()
#        dofmap = Vh.dofmap()
#        bbt = mesh.bounding_box_tree()
#        sdim = dolfin_element.space_dimension()
#        v = np.zeros(sdim)
#        
#        A = PETSc.Mat()
#        A.create(mesh.mpi_comm())
#        A.setSizes([ntargets, Vh.dim()])
#        A.setType("aij")
#        A.setPreallocationNNZ(sdim*ntargets)
#        
#        # In Parallel we will need to fix the rowLGmap so that only the points in the
#        # local mesh are kept.    
#        rowLGmap = PETSc.LGMap().create(range(0,ntargets), comm = mesh.mpi_comm() )
#        colLGmap = PETSc.LGMap().create(dofmap.dofs(), comm = mesh.mpi_comm() )
#        A.setLGMap(rowLGmap, colLGmap)
#        
#        for k in range(ntargets):
#            t = targets[k,:]
#            p = dl.Point(t)
#            cell_id = bbt.compute_first_entity_collision(p)
#            tvert = coords[cells[cell_id,:],:]
#            dolfin_element.evaluate_basis_all(v,t,tvert, cell_id)
#            cols = dofmap.cell_dofs(cell_id)
#            for j in range(sdim):
#                A[k, cols[j]] = v[j]
#                
#        A.assemblyBegin()
#        A.assemblyEnd()
#         
#        return dl.Matrix( dl.PETScMatrix(A) )