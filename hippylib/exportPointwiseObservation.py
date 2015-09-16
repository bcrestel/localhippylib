import numpy as np

def exportPointwiseObservation(points, data, fname, varname="observation"):
    """
    This function write a VTK PolyData file to visualize pointwise data.
    Inputs:
    - points:  locations of the points 
               (numpy array of size number of points by space dimension)
    - data:    pointwise values
               (dolfin vector of size number of points)
    - fname:   filename for the vtk file to export
    - varname: name of the variable for the vtk file
    """
    ndim = points.shape[1]
    npoints = points.shape[0]
    
    if ndim == 2:
        points3d = np.zeros((npoints,3), dtype = points.dtype)
        points3d[:,:-1] = points
        exportPointwiseObservation(points3d, data, fname, varname)
        return
    
    f = open(fname, 'w')
    f.write('<VTKFile version="0.1" byte_order="LittleEndian" type="PolyData">\n')
    f.write('<PolyData>\n')
    f.write('<Piece NumberOfPoints="{0:d}" NumberOfVerts="{1:d}">\n'.format(npoints,npoints))
    f.write('<Points>\n')
    f.write('<DataArray NumberOfComponents="{0:d}" format="ascii" type="Float32">\n'.format(ndim))
    f.write( np.array_str(points).replace("[", "").replace("]", "") )
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
    f.write(np.array_str( data.array() ).replace("[", "").replace("]", "") )
    f.write('\n</DataArray>\n</PointData>\n<CellData>\n</CellData>\n</Piece>\n</PolyData>\n</VTKFile>')
