"""
Script that reads results folders and creates comparisons between all results.

This script will generate:
1. A table containing the mean differences of all pairs of experiments
3. A .vtk for each pair of experiment containing the similarities at each vertex.


It also looks for the ParaView API and make calls to create the
vectorfield for each case, with the
adecuate scale and everything.
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

# Path of the experiments to compare
exp_path_1 = sys.argv[1]
exp_path_2 = sys.argv[2]

exp_name_1 = os.path.basename(os.path.dirname(exp_path_1))
exp_name_2 = os.path.basename(os.path.dirname(exp_path_2))

# output folder. Will be created if it does not exist
out_path = sys.argv[3]

if not os.path.exists(out_path):
    os.makedirs(out_path)

# Dataset on guardar els resultats
# Columnes (name)

# Columnes (lists of lists)
# One for left, one for right
files = glob.glob(exp_path_1 + 'lefthippo/' + '*_values.csv') + glob.glob(exp_path_2 + 'righthippo/' + '*_values.csv')
data_r = np.zeros([len(files), len(files)])
data_l = np.zeros([len(files), len(files)])

# For each hemisphere, separately
# Base
# Load all files
# For each file
#TODO maybe parallelize this shit
# TODO include the interactions
for d in ['lefthippo', 'righthippo']:
    files = glob.glob(exp_path_1 + d + '/*_values.csv') + glob.glob(exp_path_2 + d + '/*_values.csv')
    cols = []
    for i in range(len(files)):
        # Get the name
        fi = files[i]
        # Get experiment name
        e = os.path.basename(os.path.dirname(os.path.dirname(fi)))
        # Get test name
        name_a = os.path.basename(fi)
        name_a = name_a.replace("_values.csv", "")
        df_a = pd.read_csv(fi, header=None, names=["T", "Angle", "x", "y", "z"])
        components_a = df_a.iloc[:, [2, 3, 4]].to_numpy()

        # Save name of experiment
        cols.append(e + '_' + '_' + name_a)
        print(e + '_' + '_' + name_a)
        # Here we will save all the results of this experiment
        for j in range(i, len(files)):
            name_b = os.path.basename(files[j])
            name_b = name_b.replace("_values.csv", "")
            df_b = pd.read_csv(files[j], header=None, names=["T", "Angle", "x", "y", "z"])
            components_b = df_b.iloc[:, [2, 3, 4]].to_numpy()

            # Compute distance and compute mean
            result_similarity = [cosine_similarity(x.reshape(1, -1), y.reshape(1, -1)) for (x, y) in zip(components_a, components_b)]

            # Get median of the cosine similarity
            mean_res = np.mean(result_similarity)
            # Save result
            if d == 'lefthippo':
                data_l[i, j] = mean_res
                data_l[j, i] = mean_res
            else:
                data_r[i, j] = mean_res
                data_r[j, i] = mean_res
            # If we want the general result, use the invididual script

# Create the dataframe from the dictionary lists
df_r = pd.DataFrame(data_r, columns=cols, index=cols)
df_l = pd.DataFrame(data_l, columns=cols, index=cols)

# Save all the files
df_r.to_csv(out_path + 'dist_experiments_r.csv', index=True)
df_l.to_csv(out_path + 'dist_experiments_l.csv', index=True)
