"""
This file contains util functions for the execution of the jupyter notebooks.

In case any of the notebooks is to be used independently outside this directory, the corresponding
functions should be copied to the notebook.
"""
import vtk
import numpy as np

def ImportVTKPoints(pts, mesh):
    """
    Iserts numpy array of points into vtkPolyData
    pts numpy.array: Nx3 array of vertices
    Return updated mesh
    """
    (nv,mv) = pts.shape
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(nv)
    for j in range(nv):
        points.SetPoint( j, pts[j,0], pts[j,1], pts[j,2] )
    mesh.SetPoints(points)
    return mesh
    
    
def ExtractVTKPoints( mesh ):
    """
    Extract points from vtk structures

    mesh vtk object: mesh

    Return the Nx3 numpy.array of the vertices.
    """
    n = mesh.GetNumberOfPoints();
    vertex = np.zeros( ( n, 3 ) )
    for i in range( n ):
        mesh.GetPoint( i, vertex[i,:] )
    
    return vertex


def LoadVTKMesh(f, flatten=True):
    """
    Load a mesh from a file.
    
    filename: str with the path to a mesh
    
    Return the Nx3 numpy array of vertices.
    """
    rd = vtk.vtkPolyDataReader()
    rd.SetFileName(f)
    rd.Update()
    pts = ExtractVTKPoints(rd.GetOutput())
    n_pts = pts.shape[0]
    if flatten:
        pts = pts.flatten()
    return pts, n_pts

def apply_scalar_and_save(base, scalar, out):
    """
    Apply scalar map to a mesh.
    base: filename of the base mesh.
    scalar: the scalar map to apply
    out: where to save the result
    """
    rd = vtk.vtkPolyDataReader()
    rd.SetFileName(base)
    rd.Update()
    mesh = rd.GetOutput()
    mesh.GetPointData().SetScalars(scalar)

    rd2 = vtk.vtkPolyDataWriter()
    rd2.NewInstance()
    rd2.SetInputData(mesh)
    rd2.SetFileName(out)
    rd2.Write()