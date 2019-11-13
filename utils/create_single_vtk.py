"""
Auxiliar function to apply scalar map to a mesh. To call
from matlab. This applies 1 scalar map

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
vtk_mesh0.SetName("Pvalues")
for x in t_values:
    vtk_mesh0.InsertNextValue(x)
   

# Save it to disk
rd = vtk.vtkPolyDataReader()
rd.SetFileName(avg_mesh)
rd.Update()
mesh = rd.GetOutput()
mesh.GetPointData().SetScalars(vtk_mesh0)

rd2 = vtk.vtkPolyDataWriter()
rd2.NewInstance()
rd2.SetInputData(mesh)
rd2.SetFileName(out_file)
rd2.Write()