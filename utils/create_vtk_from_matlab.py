"""
Auxiliar function to apply scalar map to a mesh. To call
from matlab. This applies 4 scalar maps:
    % 1: T-statistic (Hotelling's T)
    % 2: T values for x coordinate
    % 3: T values for y coordinate
    % 4: T values for z coordinate

avg_mesh: path to avg mesh (string)
csv_file: path to csv file with the t values (string)
out_file: output path for the mesh (string)
"""

import sys
import vtk
import numpy as np
import pandas as pd

avg_mesh = sys.argv[1]
csv_file = sys.argv[2]
out_file = sys.argv[3]

# Hotelling T
df_values = pd.read_csv(csv_file, header=None)
t_values = df_values[0].values
vtk_mesh0 = vtk.vtkDoubleArray()
vtk_mesh0.SetName("Hotelling's T")
for x in t_values:
    vtk_mesh0.InsertNextValue(x)

t_values = df_values[1].values
vtk_mesh1 = vtk.vtkDoubleArray()
vtk_mesh1.SetName("angle")
for x in t_values:
    vtk_mesh1.InsertNextValue(x)

t_values = df_values[2].values
vtk_mesh2 = vtk.vtkDoubleArray()
vtk_mesh2.SetName("T-x")
for x in t_values:
    vtk_mesh2.InsertNextValue(x)

t_values = df_values[3].values
vtk_mesh3 = vtk.vtkDoubleArray()
vtk_mesh3.SetName("T-y")
for x in t_values:
    vtk_mesh3.InsertNextValue(x)

t_values = df_values[4].values
vtk_mesh4 = vtk.vtkDoubleArray()
vtk_mesh4.SetName("T-z")
for x in t_values:
    vtk_mesh4.InsertNextValue(x)

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

rd2 = vtk.vtkPolyDataWriter()
rd2.NewInstance()
rd2.SetInputData(mesh)
rd2.SetFileName(out_file)
rd2.Write()
