"""
Compare similarity between two meshes.

Given two meshes, with information of T-values and effects in three dimensions, 
compute a distance measure between them.
Example:
python spatial_distance -AD+CN_values.csv 
-AD+CN_values.csv test.vtk 
out_mesh_mean_l.vtk
"""

import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import sys
import vtk

# Path of mesh 1
mesh_1 = sys.argv[1]
# Path of mesh 2
mesh_2 = sys.argv[2]

out_file = sys.argv[3]

avg_mesh = sys.argv[4]

# Load the two meshes
df_values_mesh_1 = pd.read_csv(mesh_1, header=None)
df_values_mesh_2 = pd.read_csv(mesh_2, header=None)

# Remember that those meshes are organized as:
# 0 Hotellings's T
# 1 angle
# 2 x component
# 3 y component
# 4 z component
components_1 = df_values_mesh_1.iloc[:, [2, 3, 4]].to_numpy()
components_2 = df_values_mesh_2.iloc[:, [2, 3, 4]].to_numpy()

# Compute cosine similarity
result_similarity = [cosine_similarity(x.reshape(1,-1), y.reshape(1,-1)) for (x, y) in zip(components_1, components_2)]

# Get median of the cosine similarity
mean_res = np.mean(result_similarity)
print('Mean similarity: ' + str(mean_res))


# Save as resulting mesh
t_values = result_similarity
vtk_mesh = vtk.vtkDoubleArray()
vtk_mesh.SetName("Similarity")
for x in t_values:
    vtk_mesh.InsertNextValue(x)

# Save it to disk
rd = vtk.vtkPolyDataReader()
rd.SetFileName(avg_mesh)
rd.Update()
mesh = rd.GetOutput()
mesh.GetPointData().SetScalars(vtk_mesh)
rd2 = vtk.vtkPolyDataWriter()
rd2.NewInstance()
rd2.SetInputData(mesh)
rd2.SetFileName(out_file)
rd2.Write()
