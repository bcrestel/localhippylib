import dolfin as dl
import numpy as np
#import scipy.sparse as scs

try:
    from petsc4py import PETSc

    def assemblePointwiseObservation(Vh, targets):
        """
        Assemble the pointwise observation matrix:
        Input
        - Vh: FEniCS finite element space
        - targets: observation points (numpy array)
         
        Note: This Function will not work in parallel!!!
        """
        ntargets, dim = targets.shape
        mesh = Vh.mesh()
        coords = mesh.coordinates()
        cells = mesh.cells()
        dolfin_element = Vh.dolfin_element()
        dofmap = Vh.dofmap()
        bbt = mesh.bounding_box_tree()
        sdim = dolfin_element.space_dimension()
        v = np.zeros(sdim)
        
        A = PETSc.Mat()
        A.create(mesh.mpi_comm())
        A.setSizes([ntargets, Vh.dim()])
        A.setType("aij")
        A.setPreallocationNNZ(sdim*ntargets)
        
        # In Parallel we will need to fix the rowLGmap so that only the points in the
        # local mesh are kept.    
        rowLGmap = PETSc.LGMap().create(range(0,ntargets), comm = mesh.mpi_comm() )
        colLGmap = PETSc.LGMap().create(dofmap.dofs(), comm = mesh.mpi_comm() )
        A.setLGMap(rowLGmap, colLGmap)
        
        for k in range(ntargets):
            t = targets[k,:]
            p = dl.Point(t)
            cell_id = bbt.compute_first_entity_collision(p)
            tvert = coords[cells[cell_id,:],:]
            dolfin_element.evaluate_basis_all(v,t,tvert, cell_id)
            cols = dofmap.cell_dofs(cell_id)
            for j in range(sdim):
                A[k, cols[j]] = v[j]
                
        A.assemblyBegin()
        A.assemblyEnd()
         
        return dl.Matrix( dl.PETScMatrix(A) )

except:
    
    import os
    abspath = os.path.dirname( os.path.abspath(__file__) )
    sdir = os.path.join(abspath,"AssemblePointwiseObservation")
    header_file = open(os.path.join(sdir,"AssemblePointwiseObservation.h"), "r")
    code = header_file.read()
    header_file.close()
    #check the dolfin version to decide which cpp to include
    if dl.dolfin_version()[2] == "4":
        cpp_sources = ["AssemblePointwiseObservation_v14.cpp"]
    else:
        cpp_sources = ["AssemblePointwiseObservation.cpp"]  
    cpp_module = dl.compile_extension_module(
    code=code, source_directory=sdir, sources=cpp_sources,
    include_dirs=[".",  sdir])
    
    def assemblePointwiseObservation(Vh, targets):
        #Ensure that PetscInitialize is called
        dummy = dl.assemble( dl.inner(dl.TrialFunction(Vh), dl.TestFunction(Vh))*dl.dx )
        #Call the cpp module to compute the pointwise observation matrix
        tmp = cpp_module.PointwiseObservation(Vh,targets.flatten())
        #return the matrix
        return tmp.GetMatrix()



# def assemblePointwiseObservation(targets, Vh):
#     """
#     Assemble the pointwise observation matrix:
#     - targets: observation points (numpy array)
#     - Vh: FEniCS finite element space
#     
#     Note: This Function will work only in serial!!!
#     """
#     ntargets, dim = targets.shape
#     mesh = Vh.mesh()
#     coords = mesh.coordinates()
#     cells = mesh.cells()
#     dolfin_element = Vh.dolfin_element()
#     dofmap = Vh.dofmap()
#     bbt = mesh.bounding_box_tree()
#     sdim = dolfin_element.space_dimension()
#     v = np.zeros(sdim)
#     B = uBLASSparseMatrix(ntargets,Vh.dim())
#     for k in range(ntargets):
#         t = targets[k,:]
#         p = Point(t)
#         cell_id = bbt.compute_first_entity_collision(p)
#         tvert = coords[cells[cell_id,:],:]
#         dolfin_element.evaluate_basis_all(v,t,tvert, cell_id)
#         # Find the dofs for the cell
#         cols = np.array( dofmap.cell_dofs(cell_id)[:], dtype = np.uintp)
#         B.setrow(k, cols, v)
# 
#         
#     B.apply("insert") 
#     
#     return Matrix(B)
# 
# def assemblePointwiseObservation(targets, Vh):
#     """
#     Assemble the pointwise observation matrix:
#     - targets: observation points (numpy array)
#     - Vh: FEniCS finite element space
#      
#     Note: This Function will work only in serial!!!
#     """
#     ntargets, dim = targets.shape
#     mesh = Vh.mesh()
#     coords = mesh.coordinates()
#     cells = mesh.cells()
#     dolfin_element = Vh.dolfin_element()
#     dofmap = Vh.dofmap()
#     bbt = mesh.bounding_box_tree()
#     sdim = dolfin_element.space_dimension()
#     v = np.zeros(sdim)
#     rows = np.zeros(ntargets*sdim, dtype='int')
#     cols = np.zeros(ntargets*sdim, dtype='int')  
#     vals = np.zeros(ntargets*sdim)
#     for k in range(ntargets):
#         t = targets[k,:]
#         p = Point(t)
#         cell_id = bbt.compute_first_entity_collision(p)
#         tvert = coords[cells[cell_id,:],:]
#         dolfin_element.evaluate_basis_all(v,t,tvert, cell_id)
#         jj = np.arange(sdim*k,sdim*(k+1))
#         rows[jj] = k
#         # Find the dofs for the cell
#         cols[jj] = dofmap.cell_dofs(cell_id)
#         vals[jj] = v
#          
#     ij = np.concatenate((np.array([rows]), np.array([cols])),axis=0)
#     B = sps.csr_matrix((vals, ij), shape=(ntargets,Vh.dim()))
#      
#     return B
