# %% markdown
# # Preprocessing ADNI
# Preprocess the meshes of ADNI, using the same code as before, but for the ADNI. We combine the scripts 0.8 study of the data and 2.0 process images.
# %%
import pandas as pd
import os
import vtk
import numpy as np
import sys

# Import all functions in utils/utils

sys.path.append('../')
from utils.utils import *

# Important paths

# Path where the meshes are located
meshes_base = ""

# out path for right and left for intermediate results
mesh_out_1 = ""
mesh_out_2 = ""

# out path for right, left and full meshes, for final results
mesh_aligned_r = ""
mesh_aligned_l = ""
mesh_aligned = ""

# ADNIMERGE ADNI .csv path
ADNIMERGE_path= ""

# Path to the template mesh
template_mesh = ""

# Csv out with final info of all subjects which mesh was preprocessed
out_csv = ""

df_adni = pd.read_csv(ADNIMERGE_path)

# Turn all other categories into MCI or AD

dx_dict = {'EMCI': 'LMCI',
           'SMC': 'CN'}

df_adni.replace({'DX_bl': dx_dict}, inplace=True)

# Categorize SITE, process as usual

df_adni["SITE"] = df_adni["SITE"].astype('category')
df_adni["SITE_cat"] = df_adni["SITE"].cat.codes


if not os.path.exists(mesh_aligned_r):
    os.makedirs(mesh_aligned_r)
if not os.path.exists(mesh_aligned_l):
    os.makedirs(mesh_aligned_l)
if not os.path.exists(mesh_out_1):
    os.makedirs(mesh_out_1)
if not os.path.exists(mesh_out_2):
    os.makedirs(mesh_out_2)

df_adni = df_adni[df_adni["VISCODE"] == 'bl']

rd = vtk.vtkPolyDataReader()
rd.SetFileName(template_mesh)
rd.Update()
template = rd.GetOutput()

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
# Load the data
for i, row in df_adni.iterrows():
    idp = row.PTID[0:3] + '_S_' + row.PTID[6:]
    file = meshes_base + idp + '_warped.vtk'
    if not os.path.exists(file):
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
        print('Error in ' + idp)
        print('Says that there are ' + str(nr) + 'cc.')
        continue

    mesh_1 = get_region(mesh, mesh_out_1 + 'out_mesh_' + idp + '_l.vtk', 0)
    mesh_2 = get_region(mesh, mesh_out_2 + 'out_mesh_' + idp + '_r.vtk', 1)
    X_left.append(mesh_1)
    X_right.append(mesh_2)

print(len(X_left))
# %% markdown
# Compute the volumes of each mesh
# %%
# Iterate over the rows of the saved data, get the volumes
Vol_l = []
Vol_r = []
df_adni = pd.DataFrame(row_ids)
print(len(X_right))
## Load the data
i = 0
id_list = []

for _, row in df_adni.iterrows():
    # right_volume
    idp = row.PTID[0:3] + '_S_' + row.PTID[6:]
    id_list.append(idp)

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

print(len(df_adni))
print(len(Vol_r))
df_adni["vol_r"] = Vol_r
df_adni["vol_l"] = Vol_l
df_adni['adni_id'] = id_list

df_adni.PTGENDER = pd.Categorical(df_adni.PTGENDER)
df_adni["gender_int"] = df_adni.PTGENDER.cat.codes + 1

n_apoe = []
c_apoe = []
for x in df_adni.APOE4.values:
    if x == 0:
        n_apoe.append(0)
        c_apoe.append('NC')
    elif x == 1:
        n_apoe.append(1)
        c_apoe.append('HE')
    elif x == 2:
        n_apoe.append(2)
        c_apoe.append('HO')
    else:
        print('error?')

df_adni["apoe_int"] = n_apoe
df_adni["apoe_cat"] = c_apoe

n_dx = []
for x in df_adni.DX_bl.values:
    if x == 'CN':
        n_dx.append(0)
    elif x == 'LMCI':
        n_dx.append(1)
    elif x == 'AD':
        n_dx.append(2)
    else:
        print(x)
        print('error?')

df_adni["dx_int"] = n_dx

# List of id's that present visual outliers after segmentation
# Fill the list with the ids
visual_outliers_list = []

print(len(df_adni))
print(len(X_right))

# Left mesh
i = 0
meshes_left = vtk.vtkMultiBlockDataGroupFilter()
points_left = []
for mesh in X_left:
    # Add to group
    vol_r = float(df_adni.iloc[i, :].vol_r)
    vol_l = float(df_adni.iloc[i, :].vol_l)
    idp = df_adni.iloc[i, :].PTID[0:3] + '_S_' + df_adni.iloc[i, :].PTID[6:]
    if vol_r >= 1500 and vol_r < 6000 and vol_l >= 1500 and vol_l < 6000 and idp not in visual_outliers_list:
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
    vol_r = float(df_adni.iloc[i, :].vol_r)
    vol_l = float(df_adni.iloc[i, :].vol_l)
    idp = df_adni.iloc[i, :].PTID[0:3] + '_S_' + df_adni.iloc[i, :].PTID[6:]
    if vol_r >= 1500 and vol_r < 6000 and vol_l >= 1500 and vol_l < 6000 and idp not in visual_outliers_list:
        meshes_right.AddInputData(mesh)
        # Compute mean points
        points = ExtractVTKPoints(mesh)
        points_right.append(points)
    i += 1

# For both, apply the two thresholds
df_adni = df_adni[(df_adni["vol_l"] >= 1500.0) & (df_adni["vol_l"] <= 6000.0)]
df_adni = df_adni[(df_adni["vol_r"] >= 1500.0) & (df_adni["vol_r"] <= 6000.0)]
df_adni = df_adni[~df_adni["adni_id"].isin(visual_outliers_list)]

df_adni.to_csv(out_csv, index=False, index_label=False)

# Left mesh
avg_l = np.mean(np.array(points_left), axis=0)
print(avg_l.shape)
points_l = vtk.vtkPoints()
points_l.SetNumberOfPoints(len(avg_l))
for j in range(len(avg_l)):
    points_l.SetPoint(j, avg_l[j, 0], avg_l[j, 1], avg_l[j, 2])

# Right mesh
avg_r = np.mean(np.array(points_right), axis=0)
print(avg_r.shape)
points_r = vtk.vtkPoints()
points_r.SetNumberOfPoints(len(avg_r))
for j in range(len(avg_r)):
    points_r.SetPoint(j, avg_r[j, 0], avg_r[j, 1], avg_r[j, 2])

# Left mesh
import pickle

# Compute alignment of left hippocampus
pcra = vtk.vtkProcrustesAlignmentFilter()
pcra.GetLandmarkTransform().SetModeToRigidBody()
pcra.SetInputConnection(meshes_left.GetOutputPort())
pcra.Update()

out_block = pcra.GetOutput()
n = out_block.GetNumberOfBlocks()


# Extract the output and save it to new folder
for i in range(n):
    mesh = out_block.GetBlock(i)
    # get id of said
    idp = df_adni.iloc[i,:].PTID[0:3] + '_S_' + df_adni.iloc[i,:].PTID[6:]
    # save to disk
    rd = vtk.vtkPolyDataWriter()
    f1 = mesh_aligned_l + str(idp) + '_l.vtk'
    rd.SetInputData(mesh)
    rd.SetFileName(f1)
    rd.Write()

    # SAVE IT ALSO AS .OBJ
    obj = vtk.vtkMNIObjectWriter()
    f1 = mesh_aligned_l + str(idp) + '_l.obj'
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

# Extract the output and save it to new folder
for i in range(n):
    mesh = out_block.GetBlock(i)
    # get id of said
    idp = df_adni.iloc[i,:].PTID[0:3] + '_S_' + df_adni.iloc[i,:].PTID[6:]
    # save to disk
    rd = vtk.vtkPolyDataWriter()
    f1 = mesh_aligned_r + str(idp) + '_r.vtk'
    rd.SetInputData(mesh)
    rd.SetFileName(f1)
    rd.Write()

    # SAVE IT ALSO AS .OBJ
    obj = vtk.vtkMNIObjectWriter()
    f1 = mesh_aligned_r + str(idp) + '_r.obj'
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
