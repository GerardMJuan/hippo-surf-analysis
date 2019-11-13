"""
Script that merges both ADNI and ALFA datasets in order to test them separately.

This scripts moves the necessary meshes to this directory and creates the .csv combining both. 
Note that this dataset needs to be before processing, and the processing needs to be done afterwards.

need to set export LC_ALL="C"

"""
import pandas as pd
import os
import locale
locale.setlocale(locale.LC_NUMERIC, 'C')
import vtk
import numpy as np
import sys

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

# Paths of the csv detailing the two main datasets
adni_path = ""
alfa_path = ""

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

df_adni = pd.read_csv(adni_path)
df_alfa = pd.read_csv(alfa_path)

df_alfa["id_alomar"] = df_alfa["id_alomar"].astype(int)

# Merge them. In order to merge them, we need to rename the columns that are important and the same
# PTID and id_alomar should be the same
# site should be 789, i dpres s'ha de categoritzar d'alguna forma
# APOE4 i APOE_e4_num han de ser el mateix, i s'ha de prperoessar a posterori
#PTGENDER  i gender han de ser el mateix i preprocessar espais
# educations_years i PTEDUCAT tb han de ser el mateix

df_adni = df_adni[df_adni["VISCODE"] == 'bl']

# First rename columns, then merge
df_adni.rename(columns={'PTID': 'PID', 'AGE': 'age', 'APOE4': 'APOE_e4_num', 'PTGENDER': 'gender', 'PTEDUCAT': 'education_years'}, inplace=True)
df_alfa.rename(columns={'id_alomar': 'PID'}, inplace=True)

# Overwrite other diagnosis
aux_dict = {'EMCI': 'LMCI'}
df_adni.replace({'DX_bl': aux_dict}, inplace=True)

aux_dict = {'Male': 'M',
            'Female': 'F'}
df_adni.replace({'gender': aux_dict}, inplace=True)

aux_dict = {' M': 'M',
            ' F': 'F'}
df_alfa.replace({'gender': aux_dict}, inplace=True)

df_alfa["SITE"] = 789
df_alfa["DX_bl"] = 'CN'

df_total = pd.concat([df_adni, df_alfa], axis=0, ignore_index=True)
# Keep only relevant columns
df_total = df_total[["PID", "gender", "age", "education_years", "APOE_e4_num", "DX_bl", "SITE"]]

# Categorize SITE, process as usual

df_total["SITE"] = df_total["SITE"].astype('category')
df_total["SITE_cat"] = df_total["SITE"].cat.codes

df_total = df_total[df_total["DX_bl"] == 'CN']

df_total.to_csv(out_csv, index=False, index_label=False)

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

for i, row in df_total.iterrows():
    id = str(row.PID)
    file = mallas_path_out + id + '_warped.vtk'
    if not os.path.exists(file):
        # print(file)
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


df_ids.gender = pd.Categorical(df_ids.gender)
df_ids["gender_int"] = df_ids.gender.cat.codes

c_apoe = []

for x in df_ids.APOE_e4_num.values:
    if x == 0:
        c_apoe.append('NC')
    elif x == 1:
        c_apoe.append('HE')
    elif x == 2:
        c_apoe.append('HO')
    else:
        print(x)
        print('error?')

df_ids["apoe_int"] = df_ids["APOE_e4_num"]
df_ids["apoe_cat"] = c_apoe

n_dx = []
for x in df_ids.DX_bl.values:
    if x == 'CN':
        n_dx.append(0)
    elif x == 'LMCI':
        n_dx.append(1)
    elif x == 'AD':
        n_dx.append(2)
    elif x == 'SMC':
        print('lmao')
        n_dx.append(3)
    else:
        print(x)
        print('error?')

df_ids["dx_int"] = n_dx

df_ids.to_csv(out_csv, index=False, index_label=False)

# # Visualization of the volume distribution
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
sns.set()

df_ids['Subjects'] = range(len(df_ids))
# Color between non-carriers, heterozigotes and homozigotes
fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(15, 5))

ax1 = sns.scatterplot(x="Subjects", y="vol_l",
                     data=df_ids, ax=ax1)

ax2 = sns.scatterplot(x="Subjects", y="vol_r",
                     data=df_ids, ax=ax2)

ax1.set_ylabel('Volume (mm3)')
ax2.set_ylabel('Volume (mm3)')

ax1.set_title('Volume of left hippocampus')
ax2.set_title('Volume of right hippocampus')

ax1.axhline(1500, color='r')
ax2.axhline(1500, color='r')

plt.tight_layout()
plt.show()

# List of id's that present visual outliers after segmentation
# Fill the list with the ids
visual_outliers_list = []


# For both, apply the two thresholds
df_apoe_f = df_ids[(df_ids["vol_l"] >= 1500.0) & (df_ids["vol_l"] <= 6000.0)]
df_apoe_f = df_apoe_f[(df_apoe_f["vol_r"] >= 1500.0) & (df_apoe_f["vol_r"] <= 6000.0)]

df_apoe_f = df_apoe_f[~df_apoe_f["PID"].isin(visual_outliers_list)]
# What is this column drop?
# df_apoe_f.drop(df_apoe_f.columns[[0]], axis=1, inplace=True)
df_apoe_f.to_csv(out_csv, index=False, index_label=False)
print('after filtering')
print(len(df_apoe_f))

# Group left and right meshes into a MultiBlockData

# Left mesh
i = 0
meshes_left = vtk.vtkMultiBlockDataGroupFilter()
points_left = []
for mesh in X_left:
    # Add to group
    vol_r = float(df_ids.iloc[i,:].vol_r)
    vol_l = float(df_ids.iloc[i,:].vol_l)
    idb = df_ids.iloc[i,:].PID
    if vol_r >= 1500 and vol_r < 5000 and vol_l >= 1500 and vol_l < 5000 and idb not in visual_outliers_list:
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
    idb = df_ids.iloc[i,:].PID
    if vol_r >= 1500 and vol_r < 5000 and vol_l >= 1500 and vol_l < 5000 and idb not in visual_outliers_list:
        meshes_right.AddInputData(mesh)
        # Compute mean points
        points = ExtractVTKPoints(mesh)
        points_right.append(points)
    i += 1

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
    idb = df_ids.iloc[i,:].PID

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
    id = df_apoe_f.iloc[i,:].PID
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
    id = df_apoe_f.iloc[i,:].PID
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