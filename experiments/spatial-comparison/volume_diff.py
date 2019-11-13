"""
Script that computes the volume change with respect to the average for each of the tests,
using the corresponding effect

"""
import glob
import numpy as np
import os
import pandas as pd
import sys
import re
from sklearn.metrics.pairwise import cosine_similarity
import vtk
from pathlib import Path
sys.path.append('../')
from utils.utils import *


# Path of the experiment to compute, will be saved here
exp_path = sys.argv[1]

# Load either adni or alfa path
if "adni" in exp_path:
    avg_l = ''
    avg_r = ''
else:
    avg_l = ''
    avg_r = ''

# Get the points of the mesh
# Save it to disk

avg_mesh_l = vtk.vtkPolyDataReader()
avg_mesh_l.SetFileName(avg_l)
avg_mesh_l.Update()
avg_mesh_l = avg_mesh_l.GetOutput()

avg_mesh_r = vtk.vtkPolyDataReader()
avg_mesh_r.SetFileName(avg_r)
avg_mesh_r.Update()
avg_mesh_r = avg_mesh_r.GetOutput()

avg_points_l = ExtractVTKPoints(avg_mesh_l)
avg_points_r = ExtractVTKPoints(avg_mesh_r)


# Columnes for the name of experiments 
cols = []

# volume
vol = []

# volume with respect to the average
vol_diff = []

# percentage vvolume
vol_percent = []

# Hipp
hippo = []

# Compute volume for left and right

Mass = vtk.vtkMassProperties()
Mass.SetInputData(avg_mesh_l)
volume_avg_l = Mass.GetVolume() 
cols.append('average_left')
vol.append(volume_avg_l)
vol_diff.append(0.0)
hippo.append('L')
vol_percent.append(0.0)

Mass = vtk.vtkMassProperties()
Mass.SetInputData(avg_mesh_r)
volume_avg_r = Mass.GetVolume() 
cols.append('average_right')
vol.append(volume_avg_r)
vol_diff.append(0.0)
vol_percent.append(0.0)
hippo.append('R')
# For each hemisphere, separately
# For each file
for d, avg_mesh, avg, vol_avg in zip(['lefthippo', 'righthippo', 'lefthippo_ageapoe', 'righthippo_ageapoe'],
                                     [avg_mesh_l, avg_mesh_r, avg_mesh_l, avg_mesh_r],
                                     [avg_points_l, avg_points_r, avg_points_l, avg_points_r],
                                     [volume_avg_l, volume_avg_r, volume_avg_l, volume_avg_r]):
    
    # Append volume of the average, to compare
    files = glob.glob(exp_path + d + '/*_values.csv')
    for f in files:

        name = os.path.basename(f)
        name = name.replace("_values.csv", "")
        new_name = name

        df_a = pd.read_csv(f, header=None, names=["T", "Angle", "x", "y", "z"])
        components_a = df_a.iloc[:, [2, 3, 4]].to_numpy()

        # Apply components 
        new_points = avg + components_a
        points_mesh = vtk.vtkPoints()
        points_mesh.SetNumberOfPoints(len(new_points))
        for j in range(len(new_points)):
            points_mesh.SetPoint(j, new_points[j, 0], new_points[j, 1], new_points[j, 2])

        # Get the mean mesh and save it
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points_mesh)
        polydata.SetPolys(avg_mesh.GetPolys())
        polydata.Modified()

        Mass = vtk.vtkMassProperties()
        Mass.SetInputData(polydata)
        volume_new = Mass.GetVolume()

        # Add to structure
        cols.append(new_name)
        vol.append(volume_new)
        vol_diff.append(volume_new - vol_avg)
        vol_percent.append(((volume_new - vol_avg) / vol_avg) * 100)
        if d == 'lefthippo':
            hippo.append('L')
        else:
            hippo.append('R')

# Create the dataframe from the dictionary lists
dict_out = {
    "experiment": cols,
    "H": hippo,
    "volume": vol,
    "Volume diff": vol_diff,
    "Vol %": vol_percent
}
df_base = pd.DataFrame(dict_out)

# Save all the files
df_base.to_csv(exp_path + 'volume_table.csv', index=False)