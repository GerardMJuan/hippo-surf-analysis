"""
Auxiliar function to apply scalar map to a mesh. This is a variance
of create_vtk_from_matlab, which 
    % 1: T values for x coordinate, linear
    % 2: T values for y coordinate, linear
    % 3: T values for z coordinate, linear
    % 4: T values for x coordinate, squared
    % 5: T values for y coordinate, squared
    % 6: T values for z coordinate, squared
avg_mesh: path to avg mesh (string)
csv_file: path to csv file with the t values (string)
out_file: output path for the mesh (string)
"""

import sys
import vtk
import pandas as pd

avg_mesh = sys.argv[1]
csv_file = sys.argv[2]
out_file = sys.argv[3]

# Hotelling T
df_values = pd.read_csv(csv_file, header=None)
t_values = df_values[0].values
vtk_mesh0 = vtk.vtkDoubleArray()
vtk_mesh0.SetName("effxlin")
for x in t_values:
    vtk_mesh0.InsertNextValue(x)

t_values = df_values[1].values
vtk_mesh1 = vtk.vtkDoubleArray()
vtk_mesh1.SetName("effylin")
for x in t_values:
    vtk_mesh1.InsertNextValue(x)

t_values = df_values[2].values
vtk_mesh2 = vtk.vtkDoubleArray()
vtk_mesh2.SetName("effzlin")
for x in t_values:
    vtk_mesh2.InsertNextValue(x)

t_values = df_values[3].values
vtk_mesh3 = vtk.vtkDoubleArray()
vtk_mesh3.SetName("effxsq")
for x in t_values:
    vtk_mesh3.InsertNextValue(x)

t_values = df_values[4].values
vtk_mesh4 = vtk.vtkDoubleArray()
vtk_mesh4.SetName("effysq")
for x in t_values:
    vtk_mesh4.InsertNextValue(x)

t_values = df_values[5].values
vtk_mesh5 = vtk.vtkDoubleArray()
vtk_mesh5.SetName("effzsq")
for x in t_values:
    vtk_mesh5.InsertNextValue(x)

# Save it to disk
rd = vtk.vtkPolyDataReader()
rd.SetFileName(avg_mesh)
rd.Update()


mesh = rd.GetOutput()
mesh.GetPointData().SetScalars(vtk_mesh0)
mesh.GetPointData().AddArray(vtk_mesh1)
mesh.GetPointData().AddArray(vtk_mesh2)
mesh.GetPointData().AddArray(vtk_mesh3)
mesh.GetPointData().AddArray(vtk_mesh4)
mesh.GetPointData().AddArray(vtk_mesh5)

rd2 = vtk.vtkPolyDataWriter()
rd2.NewInstance()
rd2.SetInputData(mesh)
rd2.SetFileName(out_file)
rd2.Write()
