# # Process images without resizing.
# date: 7-Jan-2019
# We will process the images with procrustes analysis, but without resizing.
# We will also try to do the analysis without dividing the hippocampus between left
# and right, to see what happens. Reusing most of the code from alicia_reproduce/0.5-Preprocessing
# %%
import pandas as pd
import os
import locale
locale.setlocale(locale.LC_NUMERIC, 'C')
import vtk
import numpy as np
import sys

locale.setlocale(locale.LC_NUMERIC, 'C')

# Import all functions in utils/utils
sys.path.append('../')
from utils.utils import *

# Path where the meshes are
mallas_path_out = ""

# out path for right and left for intermediate results
mesh_out_1 = ""
mesh_out_2 = ""

# out path for right, left and full meshes, for final results
mesh_aligned_r = ""
mesh_aligned_l = ""
mesh_aligned = ""

# .ALFA csv info path
ALFA_path = ""

# Path to the template mesh
template_mesh = ""

# Csv out with final info of all subjects which mesh was preprocessed
out_csv = ""

if not os.path.exists(mesh_aligned_r):
    os.makedirs(mesh_aligned_r)
if not os.path.exists(mesh_aligned_l):
    os.makedirs(mesh_aligned_l)
if not os.path.exists(mesh_aligned):
    os.makedirs(mesh_aligned)

if not os.path.exists(mesh_out_1):
    os.makedirs(mesh_out_1)
if not os.path.exists(mesh_out_2):
    os.makedirs(mesh_out_2)


df_apoe = pd.read_csv(ALFA_path)

print(len(df_apoe))
df_apoe = df_apoe.dropna()
print(len(df_apoe))
df_apoe.head()


rd = vtk.vtkPolyDataReader()
rd.SetFileName(template_mesh)
rd.Update()
template = rd.GetOutput()
# %% markdown
# First, we need to divide each mesh between left and right. To do this, we use a VTK function to divide each mesh into connected components.
# %%
# auxiliary function to not repeat code
def get_region(mesh, mesh_out, region):
    """
    From a region with two connected components,
    get either first or second connected component
    and save them to mesh_out.

    Returns the output mesh
    """
    # Get region
    cc = vtk.vtkPolyDataConnectivityFilter()
    cc.SetInputData(mesh)
    cc.SetExtractionModeToSpecifiedRegions()
    cc.AddSpecifiedRegion(region)
    cc.Update()
    mesh_1 = cc.GetOutput(0)

    # Remove connected points
    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(mesh_1)
    clean.Update()
    mesh_1 = clean.GetOutput(0)

    # Save
    rd = vtk.vtkPolyDataWriter()
    rd.SetInputData(mesh_1)
    rd.SetFileName(mesh_out)
    rd.Write()

    return mesh_1

# We store the data here
X_right = []
X_left = []
X_total = []

# actual ids
row_ids = []

# Number of points per mesh
n_pts = 0

## Load the data
for i, row in df_apoe.iterrows():
    id = str(int(row.id_alomar))
    file = mallas_path_out + id + '_warped.vtk'
    if not os.path.exists(file):
        print(file)
        continue

    # Add row to df_ids in order
    row_ids.append(row)
    rd = vtk.vtkPolyDataReader()
    rd.SetFileName(file)
    rd.Update()
    mesh = rd.GetOutput()
    # Add the full mesh
    X_total.append(mesh)
    cc = vtk.vtkPolyDataConnectivityFilter()
    cc.SetInputData(mesh)
    cc.SetExtractionModeToAllRegions()
    cc.Update()
    nr = cc.GetNumberOfExtractedRegions()
    # Confirm there are only 2 regions
    if nr != 2:
        print('Error in ' + id)
        print('Says that there are ' + str(nr) + 'cc.')
        continue

    mesh_1 = get_region(mesh, mesh_out_1 + 'out_mesh_' + id + '_l.vtk', 0)
    mesh_2 = get_region(mesh, mesh_out_2 + 'out_mesh_' + id + '_r.vtk', 1)
    X_left.append(mesh_1)
    X_right.append(mesh_2)

print(len(X_left))
# %% markdown
# Then, compute the volumes of each connected component, compare the volume to the volumes registered on the .csv and validate if they are the same
# %%
# Iterate over the rows of the saved data, get the volumes
Vol_l = []
Vol_r = []
df_ids = pd.DataFrame(row_ids)
print(len(df_ids))
print(len(X_right))
## Load the data
i = 0
for _, row in df_ids.iterrows():
    # right_volume
    Mass = vtk.vtkMassProperties()
    Mass.SetInputData(X_right[i])
    Mass.Update()
    Vol_r.append(Mass.GetVolume())

    # left_volume
    Mass = vtk.vtkMassProperties()
    Mass.SetInputData(X_left[i])
    Mass.Update()
    Vol_l.append(Mass.GetVolume())
    i = i + 1

df_ids["vol_r"] = Vol_r
df_ids["vol_l"] = Vol_l

df_ids.to_csv(out_csv)

df_apoe = pd.read_csv(out_csv)
df_apoe = df_apoe.dropna()
df_apoe.head()

df_apoe.gender = pd.Categorical(df_apoe.gender)
df_apoe["gender_int"] = df_apoe.gender.cat.codes + 1

n_apoe = []
c_apoe = []
for x in df_apoe.apoe.values:
    if x in [22,23,33]:
        n_apoe.append(0)
        c_apoe.append('NC')
    elif x in [24,34]:
        n_apoe.append(1)
        c_apoe.append('HE')
    elif x == 44:
        n_apoe.append(2)
        c_apoe.append('HO')
    else:
        print('error?')

df_apoe["apoe_int"] = n_apoe
df_apoe["apoe_cat"] = c_apoe

# # Visualization of the volume distribution
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
sns.set()

df_apoe['Subjects'] = range(len(df_apoe))
# Color between non-carriers, heterozigotes and homozigotes
fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(15, 5))

ax1 = sns.scatterplot(x="Subjects", y="vol_l",
                     data=df_apoe, ax=ax1)

ax2 = sns.scatterplot(x="Subjects", y="vol_r",
                     data=df_apoe, ax=ax2)

ax1.set_ylabel('Volume (mm3)')
ax2.set_ylabel('Volume (mm3)')

ax1.set_title('Volume of left hippocampus')
ax2.set_title('Volume of right hippocampus')

ax1.axhline(1500, color='r')
ax2.axhline(1500, color='r')

plt.tight_layout()

# List of id's that present visual outliers after segmentation
# Fill the list with the ids
visual_outliers_list = []

# visual_outliers_list = []

# For both, apply the two thresholds
df_apoe_f = df_apoe[df_apoe["vol_l"] >= 1500.0]
df_apoe_f = df_apoe_f[df_apoe_f["vol_r"] >= 1500.0]
df_apoe_f = df_apoe_f[~df_apoe_f["id_alomar"].isin(visual_outliers_list)]
df_apoe_f.drop(df_apoe_f.columns[[0]], axis=1, inplace=True)
df_apoe_f.to_csv(out_csv, index=False, index_label=False)
print(len(df_apoe_f))
# %% markdown
# For both left and right hippocampus, use similarity alignment to alineate the meshes, using vtk.vtkProcustesAlignmentFilter()
# %%
# Group left and right meshes into a MultiBlockData
# la patillada del if vol xd
# Prova random: calcular el mean a mà, punt a punt.

# Left mesh
i = 0
meshes_left = vtk.vtkMultiBlockDataGroupFilter()
points_left = []
for mesh in X_left:
    # Add to group
    vol_r = float(df_ids.iloc[i,:].vol_r)
    vol_l = float(df_ids.iloc[i,:].vol_l)
    idb = int(df_ids.iloc[i,:].id_alomar)
    if vol_r >= 2000 and vol_l >= 1500 and idb not in visual_outliers_list:
        meshes_left.AddInputData(mesh)
        # Compute mean points
        points = ExtractVTKPoints(mesh)
        points_left.append(points)
    i += 1

# Right mesh
i = 0
meshes_right = vtk.vtkMultiBlockDataGroupFilter()
points_right = []
for mesh in X_right:
    # Add to group
    vol_r = float(df_ids.iloc[i,:].vol_r)
    vol_l = float(df_ids.iloc[i,:].vol_l)
    idb = int(df_ids.iloc[i,:].id_alomar)

    if vol_r >= 2000 and vol_l >= 1500 and idb not in visual_outliers_list:
        meshes_right.AddInputData(mesh)
        # Compute mean points
        points = ExtractVTKPoints(mesh)
        points_right.append(points)
    i += 1

# %%
# Left mesh
avg_l = np.mean(np.array(points_left), axis = 0)
print(avg_l.shape)
points_l = vtk.vtkPoints()
points_l.SetNumberOfPoints(len(avg_l))
for j in range(len(avg_l)):
    points_l.SetPoint( j, avg_l[j,0], avg_l[j,1], avg_l[j,2] )

# Right mesh
avg_r = np.mean(np.array(points_right), axis = 0)
print(avg_r.shape)
points_r = vtk.vtkPoints()
points_r.SetNumberOfPoints(len(avg_r))
for j in range(len(avg_r)):
    points_r.SetPoint( j, avg_r[j,0], avg_r[j,1], avg_r[j,2] )


# Total mesh
i = 0
meshes_total = vtk.vtkMultiBlockDataGroupFilter()
points_total = []
for mesh in X_total:
    # Add to group
    vol_r = float(df_ids.iloc[i,:].vol_r)
    vol_l = float(df_ids.iloc[i,:].vol_l)
    idb = int(df_ids.iloc[i,:].id_alomar)

    if vol_r >= 1500 and vol_l >= 1500 and idb not in visual_outliers_list:
        meshes_total.AddInputData(mesh)
        # Compute mean points
        points = ExtractVTKPoints(mesh)
        points_total.append(points)
    i += 1

# Left mesh
import pickle

# Compute alignment of left hippocampus
pcra = vtk.vtkProcrustesAlignmentFilter()
pcra.GetLandmarkTransform().SetModeToRigidBody()
pcra.SetInputConnection(meshes_left.GetOutputPort())
pcra.Update()

out_block = pcra.GetOutput()
n = out_block.GetNumberOfBlocks()
print(len(df_apoe_f))
print(n)

# Extract the output and save it to new folder
for i in range(n):
    mesh = out_block.GetBlock(i)
    # get id of said
    id = str(int(df_apoe_f.iloc[i,:].id_alomar))
    # save to disk
    rd = vtk.vtkPolyDataWriter()
    f1 = mesh_aligned_l + 'out_mesh_' + str(id) + '_l.vtk'
    rd.SetInputData(mesh)
    rd.SetFileName(f1)
    rd.Write()

    # SAVE IT ALSO AS .OBJ
    obj = vtk.vtkMNIObjectWriter()
    f1 = mesh_aligned_l + 'out_mesh_' + str(id) + '_l.obj'
    obj.SetInputData(mesh)
    obj.SetFileName(f1)
    obj.Write()

# Get the mean mesh and save it
polydata = vtk.vtkPolyData()
mean_pts_l = pcra.GetMeanPoints()

test = ExtractVTKPoints(mean_pts_l)
np.save("mean_pts_l.npy", test)

polydata.SetPoints(mean_pts_l)
# polydata.SetPoints(points_l)
polydata.SetPolys(X_left[2].GetPolys())
polydata.Modified()

f1 = mesh_aligned_l + 'out_mesh_mean_l.vtk'
rd1 = vtk.vtkPolyDataWriter()
rd1.NewInstance()
rd1.SetInputData(polydata)
rd1.SetFileName(f1)
rd1.Write()

# SAVE IT ALSO AS .OBJ
obj = vtk.vtkMNIObjectWriter()
f1 = mesh_aligned_l + 'out_mesh_mean_l.obj'
obj.NewInstance()
obj.SetInputData(polydata)
obj.SetFileName(f1)
obj.Write()
# %%
# Right mesh
import pickle

# Compute alignment of right hippocampus
pcra2 = vtk.vtkProcrustesAlignmentFilter()
pcra2.GetLandmarkTransform().SetModeToRigidBody()
pcra2.SetInputConnection(meshes_right.GetOutputPort())
pcra2.Update()

out_block = pcra2.GetOutput(0)
n = out_block.GetNumberOfBlocks()
print(n)
print(len(df_apoe_f))

# Extract the output and save it to new folder
for i in range(n):
    mesh = out_block.GetBlock(i)
    # get id of said
    id = str(int(df_apoe_f.iloc[i,:].id_alomar))
    # save to disk
    rd = vtk.vtkPolyDataWriter()
    f1 = mesh_aligned_r + 'out_mesh_' + str(id) + '_r.vtk'
    rd.SetInputData(mesh)
    rd.SetFileName(f1)
    rd.Write()

    # SAVE IT ALSO AS .OBJ
    obj = vtk.vtkMNIObjectWriter()
    f1 = mesh_aligned_r + 'out_mesh_' + str(id) + '_r.obj'
    obj.SetInputData(mesh)
    obj.SetFileName(f1)
    obj.Write()

# Get the mean mesh and save it
polydata = vtk.vtkPolyData()
mean_pts_r = pcra2.GetMeanPoints()
test = ExtractVTKPoints(mean_pts_r)
np.save("mean_pts_r.npy", test)

polydata.SetPoints(mean_pts_r)
# polydata.SetPoints(points_r)
polydata.SetPolys(X_right[2].GetPolys())
polydata.Modified()

f1 = mesh_aligned_r + 'out_mesh_mean_r.vtk'
rd2 = vtk.vtkPolyDataWriter()
rd2.NewInstance()
rd2.SetInputData(polydata)
rd2.SetFileName(f1)
rd2.Write()

# SAVE IT ALSO AS .OBJ
obj = vtk.vtkMNIObjectWriter()
f1 = mesh_aligned_r + 'out_mesh_mean_r.obj'
obj.NewInstance()
obj.SetInputData(polydata)
obj.SetFileName(f1)
obj.Write()
# %%


# Full mesh, containing both hippocampus

# Compute alignment of right hippocampus
pcra3 = vtk.vtkProcrustesAlignmentFilter()
pcra3.GetLandmarkTransform().SetModeToRigidBody()
pcra3.SetInputConnection(meshes_total.GetOutputPort())
pcra3.Update()

out_block = pcra3.GetOutput(0)
n = out_block.GetNumberOfBlocks()
print(n)
print(len(df_apoe_f))

# Extract the output and save it to new folder
for i in range(n):
    mesh = out_block.GetBlock(i)
    # get id of said
    id = str(int(df_apoe_f.iloc[i, :].id_alomar))
    # save to disk
    rd = vtk.vtkPolyDataWriter()
    f1 = mesh_aligned + 'out_mesh_' + str(id) + '_t.vtk'
    rd.SetInputData(mesh)
    rd.SetFileName(f1)
    rd.Write()

    # SAVE IT ALSO AS .OBJ
    obj = vtk.vtkMNIObjectWriter()
    f1 = mesh_aligned + 'out_mesh_' + str(id) + '_t.obj'
    obj.SetInputData(mesh)
    obj.SetFileName(f1)
    obj.Write()

# Get the mean mesh and save it
polydata = vtk.vtkPolyData()
mean_pts_t = pcra3.GetMeanPoints()
test = ExtractVTKPoints(mean_pts_t)
np.save("mean_pts_t.npy", test)

polydata.SetPoints(mean_pts_t)
polydata.SetPolys(X_total[2].GetPolys())
polydata.Modified()

f1 = mesh_aligned + 'out_mesh_mean_t.vtk'
rd2 = vtk.vtkPolyDataWriter()
rd2.NewInstance()
rd2.SetInputData(polydata)
rd2.SetFileName(f1)
rd2.Write()

# SAVE IT ALSO AS .OBJ
obj = vtk.vtkMNIObjectWriter()
f1 = mesh_aligned + 'out_mesh_mean_t.obj'
obj.NewInstance()
obj.SetInputData(polydata)
obj.SetFileName(f1)
obj.Write()
# %%
