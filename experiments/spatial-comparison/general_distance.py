"""
General script that perform all distance comparisons that are needed.
We basically compare.
1. ALFA baseline effects with ADNI baseline effects.
2. ALFA interaction effects apoe x age with ADNI baseline effects
3. ALFA interaction effects apoe x age with ADNI interaction effects with age and apoe
4. ALFA baseline effects with ADNI NC baseline effects.
5. ALFA interaction effects with ADNI NC interaction effects.
6. ALFA baseline effects with ALL baseline effects.
6. ALFA interaction effects with ALL interaction effects.

All comparison are appropiatevly saved in a defined directory

"""

# imports
from spatial_distance import test_cosine_distance
import numpy as np
import pandas as pd
import os
import sys
import matlab.engine

#Auxiliar function

def fwer_correction(work_dir, csv_file, avg_file, save_file):
    """
    Auxiliar function to run matlab code.

    Call matlab function
    """
    eng = matlab.engine.start_matlab()
    eng.addpath('surfstat/')
    eng.addpath('experiments/spatial-comparison/')
    eng.correct_pval(work_dir, csv_file, avg_file, save_file, nargout=0)


## Auxiliar dictionaries
# pvpython = "C:/Users/gerar/Documents/CODE/pvpython"
pvpython = '/home/gerard/Documents/LIB/ParaView/bin/pvpython'

# Dict for left and right hippo:
d_hip = {
        "lefthippo": "L",
        "righthippo": "R",
        "lefthippo_ageapoe": "L",
        "righthippo_ageapoe": "R"
        }

# Dict for the test names
# will use regular expressions
d_test_names = {
    "+age": "Age",
    "-age": "Age",
    "yed": "Years-of-education",
    "-yed": "Years-of-education",
    "+volume": "Volume",
    "-volume": "Volume",
    "site": "Site",
    "-site": "Site",
    "+male_-female": "Sex",
    "+female_-male": "Sex",
    "+HO_-HE": "HO>HE",
    "+HE_-HO": "HO>HE",
    "HE+NC-HO": "APOE-Recessive",
    "-HE-NC+HO": "APOE-Recessive",
    "+NC_-HO": "HO>NC",
    "+HO_-NC": "HO>NC",
    "-NC_HE_HO": "APOE-Additive",
    "NC_HE_HO": "APOE-Additive",
    "-HE-HO+NC": "APOE-Dominant",
    "HE+HO-NC": "APOE-Dominant",
    "+HE_-NC": "HE>NC",
    "+NC_-HE": "HE>NC",
    "-AD+MCI": "AD>MCI",
    "+AD-MCI": "AD>MCI",
    "+CN-LMCI": "MCI>CN",
    "+MCI-CN": "MCI>CN",
    "CN+MCI-AD": "CN+MCI>AD",
    "CN-MCI+AD": "AD+MCI>CN",
    "AD_MCI_CN": "AD_MCI_CN",
    "-AD_MCI_CN": "AD_MCI_CN",
    "+AD-CN": "AD>CN",
    "-AD+CN": "AD>CN",
    "+AD_-CN": "AD>CN",
    "+AD_-MCI": "AD>MCI",
    "+CN_-AD": "AD>CN",
    "+CN_-MCI": "MCI>CN",
    "+MCI_-AD": "AD>MCI",
    "+MCI_-CN": "MCI>CN",
    "-CN_MCI_AD": "AD_MCI_CN",
    "CN_MCI_AD": "AD_MCI_CN",
    "-MCI-CN+AD": "CN+MCI>AD",
    "-MCIAD_CN": "AD+MCI>CN",
    "MCI+CN-AD": "CN+MCI>AD",
    "MCIAD_CN": "AD+MCI>CN",
}

# for the interaction
d_test_names_i = {
    "agesq": "Age^2x",
    "age": "Agex"
}

#Dictionary for the paraview files
paraview_dict = {"L": ("/home/gerard/Documents/CODE/nl-hippo-apoe/figures/Paraview_fig_pvalues/paraview_left_pvalue_1.py",
                               "/home/gerard/Documents/CODE/nl-hippo-apoe/figures/Paraview_fig_pvalues/paraview_left_pvalue_2.py"),
                 "R": ("/home/gerard/Documents/CODE/nl-hippo-apoe/figures/Paraview_fig_pvalues/paraview_right_pvalue_1.py",
                                "/home/gerard/Documents/CODE/nl-hippo-apoe/figures/Paraview_fig_pvalues/paraview_right_pvalue_2.py")}

#Dictionary for the paraview files
paraview_dict_sim = {"L": ("/home/gerard/Documents/CODE/nl-hippo-apoe/figures/Paraview_fig_sim/left_1.py",
                               "/home/gerard/Documents/CODE/nl-hippo-apoe/figures/Paraview_fig_sim/left_2.py"),
                 "R": ("/home/gerard/Documents/CODE/nl-hippo-apoe/figures/Paraview_fig_sim/right_1.py",
                                "/home/gerard/Documents/CODE/nl-hippo-apoe/figures/Paraview_fig_sim/right_2.py")}


#dictionary for the average mean 
avg_dict = {"L": "/home/gerard/Documents/DATA/Hippo_ALLDATA/Mallas/out_noresize_l/out_mesh_mean_l.vtk",
            "R": "/home/gerard/Documents/DATA/Hippo_ALLDATA/Mallas/out_noresize_r/out_mesh_mean_r.vtk"}

#dictionary for the average mean (obj)
avg_dict_obj = {"L": "/home/gerard/Documents/DATA/Hippo_ALLDATA/Mallas/out_noresize_l/out_mesh_mean_l.obj",
                "R": "/home/gerard/Documents/DATA/Hippo_ALLDATA/Mallas/out_noresize_r/out_mesh_mean_r.obj"}

# Change directory 
base_directory = "/home/gerard/Documents/EXPERIMENTS/hippomultivariate-FINAL-yed"
output_directory = "/home/gerard/Documents/EXPERIMENTS/hippomultivariate-FINAL-yed-sim"
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

def turn_test_name(test_name):
    # get test name
    # We are going to do a try/except: if the full string is in the dictionary, cool. 
    # If not, iterate
    try:
        test_name = d_test_names[test_name]
    except Exception:
        # Change the name to the actua ldictionary
        for item in d_test_names.keys():
            if item in test_name:
                test_name = test_name.replace(item, d_test_names[item])

        # This will only happen if the previous fails
        for item in d_test_names_i.keys():
            if item in test_name:
                test_name = test_name.replace(item, d_test_names_i[item])
    return test_name

# Tests function
def run_tests(out_dir, side_folders, dir_name_1, dir_name_2, list_tests):
    """
    Run the set of tests.
    out_dir is the final directory where all will be stored
    side_folders are the side folder (left/right). Two tuples.
    dir_name_1 and dir_name_2 are the folders where the experiments are
    list_tests are the lists of tests to do
    """
    dict_clusters = []
    dict_base = []

    for test in list_tests:
        #get new names
        test_new = (turn_test_name(test[0]), turn_test_name(test[1]))
        for d in side_folders:
            print(d)
            # Get the two files for that experiments
            test1 = dir_1 + d[0] + '/' + test[0] + '_values.csv'
            test2 = dir_2 + d[1] + '/' + test[1] + '_values.csv'
            print(test)
            # define output file 
            out_file = out_dir + d_hip[d[0]] + '-' + test[0] + '-' + test[1]

            #Run the test
            
            ## ALGUNA COSA FALLA A L'HORA DE GUARDAR LES COSES AMB ELS SIMBOLS _ < <_
            mean_res, out_csv_file = test_cosine_distance(test1, test2, out_file, avg_dict[d_hip[d[0]]], True)
            print("FWER correction")
            # here we would run the MATLAB correction
            fwer_correction(out_dir, os.path.basename(out_csv_file), avg_dict_obj[d_hip[d[0]]], os.path.basename(out_file) + '_corr')

            # here we would call paraview for automatic visualization
            fig = out_file + '_pval.vtk'
            out_fig_1 = out_file + '_pvfields_1.png'
            out_fig_2 = out_file + '_pvfields_2.png'
            print("Run paraview.")
            # Run both
            cmd = [pvpython, paraview_dict[d_hip[d[0]]][0], fig, out_fig_1]
            print(" ".join(cmd))
            os.system(" ".join(cmd))
            # Run both
            print(" ".join(cmd))
            cmd = [pvpython, paraview_dict[d_hip[d[0]]][1], fig, out_fig_2]
            os.system(" ".join(cmd))

            # here we would call paraview for automatic visualization
            fig = out_file + '.vtk'
            out_fig_1 = out_file + '_sim_1.png'
            out_fig_2 = out_file + '_sim_2.png'
            print("Run paraview sim.")
            # Run both
            cmd = [pvpython, paraview_dict_sim[d_hip[d[0]]][0], fig, out_fig_1]
            print(" ".join(cmd))
            os.system(" ".join(cmd))
            # Run both
            print(" ".join(cmd))
            cmd = [pvpython, paraview_dict_sim[d_hip[d[0]]][1], fig, out_fig_2]
            os.system(" ".join(cmd))

            # here we would save all the information to the table 

            # Open the corresponding clusters and peaks (if applicable) and save that information too.
            try:
                print(out_file + '_clusters.txt')
                print(out_file + '_peak.txt')
                df_clus = pd.read_csv(out_file + '_corr_clusters.txt')
                df_peaks = pd.read_csv(out_file + '_corr_peak.txt')
                Nclus = len(df_clus)
                if len(df_clus) > 0:
                    for row_clus in df_clus.itertuples():
                        dict_clusters.append({
                            "Test": test_new[0] + "-" + test_new[1],
                            "Hipp": d_hip[d[0]],
                            "Clus": row_clus.clusid,
                            "N": row_clus.nverts,
                            "P_cluster": row_clus.pmean
                        })
                        print(dict_clusters)
            except FileNotFoundError:
                Nclus = 0

            dict_base.append({
                "Test": test_new[0] + "-" + test_new[1],
                "Hipp": d_hip[d[0]],
                "avg_mean": mean_res,
                "Nclusters": Nclus
            })

    # columns are: test (name of test) Hipp (L or R), avg-T2 (over the surface) max-T2, and nclusters
    columns = ["Test", "Hipp", "avg_mean", "Nclusters"]

    # Columns for cluster
    columns_clus = ["Test", "Hipp", "Clus", "N", "P_cluster"]

    ## Save to disk, in both latex and csv form
    # Create the dataframe from the dictionary lists
    df_base = pd.DataFrame(dict_base, columns=columns)
    df_base.sort_values(by=["Test", "Hipp"], ascending=[True, False], inplace=True)
    df_base.drop_duplicates(inplace=True)

    df_clusters = pd.DataFrame(dict_clusters, columns=columns_clus)
    df_clusters.sort_values(by=["Test", "Hipp"], ascending=[True, False], inplace=True)
    df_clusters.drop_duplicates(inplace=True)

    # Save all the files

    df_base.to_csv(out_dir + 'base_table.csv', index=False)
    df_base.to_latex(out_dir + 'base_table.tex')

    df_clusters.to_csv(out_dir + 'clusters_table.csv', index=False)
    df_clusters.to_latex(out_dir + 'clusters_table.tex')

#####################################################
### BLOCKS OF COMPARISONS/TESTS START HERE ##########
#####################################################


"""
#####################################################
#0. ALFA baseline effects with ALFA baseline effects.
#####################################################
dir_1 = base_directory + "/alfa_baseline_yed_novol/"
dir_2 = base_directory + "/alfa_baseline_yed_novol/"
out_dir = output_directory + "/0_ALFA_ALFA_baseline/"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Define information on the directories and the tests to do
dir_name_1 = "ALFA"
dir_name_2 = "ALFA"
alfa_adni_baseline_tests = [
    ("+age", "+age"),
    ("+female_-male", "+female_-male"),
    ("NC_HE_HO", "NC_HE_HO"),
    ("HE+HO-NC", "HE+HO-NC"),
    ("NC_HE_HO", "AD_MCI_CN"),
    ("NC_HE_HO", "+age"),
    ("NC_HE_HO", "-age"),
    ("NC_HE_HO", "yed")]


# Folders of the hippocampal sides
hippo_folders = [('lefthippo', 'lefthippo'),('righthippo', 'righthippo')]


# Run the actual tests
run_tests(out_dir, hippo_folders, dir_name_1, dir_name_2, alfa_adni_baseline_tests)
"""

#####################################################
#1. ALFA baseline effects with ADNI baseline effects.
#####################################################
dir_1 = base_directory + "/alfa_baseline_yed_novol/"
dir_2 = base_directory + "/adni_baseline_yed_novol/"
out_dir = output_directory + "/1_ALFA_ADNI_baseline/"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Define information on the directories and the tests to do
dir_name_1 = "ALFA"
dir_name_2 = "ADNI"
alfa_adni_baseline_tests = [
    ("+age", "+age"),
    ("+female_-male", "+female_-male"),
    ("NC_HE_HO", "NC_HE_HO"),
    ("HE+HO-NC", "HE+HO-NC"),
    ("NC_HE_HO", "AD_MCI_CN"),
    ("NC_HE_HO", "+age"),
    ("NC_HE_HO", "-age"),
    ("HE+HO-NC", "+age"),
    ("NC_HE_HO", "yed")]


# Folders of the hippocampal sides
hippo_folders = [('lefthippo', 'lefthippo'),('righthippo', 'righthippo')]


# Run the actual tests
run_tests(out_dir, hippo_folders, dir_name_1, dir_name_2, alfa_adni_baseline_tests)

########################################################################
#2. ALFA interaction lin effects apoe x age with ADNI baseline effects.#
########################################################################
dir_1 = base_directory + "/alfa_interaction_test_lin_novol_yed/"
dir_2 = base_directory + "/adni_baseline_yed_novol/"
out_dir = output_directory + "/2_ALFA_int_lin_ADNI_baseline/"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

dir_name_1 = "ALFA-int"
dir_name_2 = "ADNI"


tests = [
    ("ageNC_HE_HO", "NC_HE_HO"),
    ("ageHE+HO-NC", "HE+HO-NC"),
    ("ageNC_HE_HO", "AD_MCI_CN"),
    ("ageNC_HE_HO", "+AD-CN"),
    ("ageNC_HE_HO", "+age"),
    ("ageNC_HE_HO", "-age"),
    ("ageHE+HO-NC", "+age")]

# Folders of the hippocampal sides
hippo_folders = [('lefthippo_ageapoe', 'lefthippo'),('righthippo_ageapoe', 'righthippo')]

# Run the actual tests
run_tests(out_dir, hippo_folders, dir_name_1, dir_name_2, tests)

########################################################################
#2.1. ALFA interaction sq effects apoe x age with ADNI baseline effects.#
########################################################################
dir_1 = base_directory + "/alfa_interaction_test_novol_yed/"
dir_2 = base_directory + "/adni_baseline_yed_novol/"
out_dir = output_directory + "/2_2_ALFA_int_ADNI_baseline/"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

dir_name_1 = "ALFA-int"
dir_name_2 = "ADNI"

tests = [
    ("agesqNC_HE_HO", "NC_HE_HO"),
    ("agesqHE+HO-NC", "HE+HO-NC"),
    ("agesqNC_HE_HO", "AD_MCI_CN"),
    ("agesqNC_HE_HO", "+AD-CN"),
    ("agesqNC_HE_HO", "+age"),
    ("agesqHE+HO-NC", "+age")]



# Folders of the hippocampal sides
hippo_folders = [('lefthippo_ageapoe', 'lefthippo'), ('righthippo_ageapoe', 'righthippo')]

# Run the actual tests
run_tests(out_dir, hippo_folders, dir_name_1, dir_name_2, tests)
"""
########################################################################################
#3. ALFA interaction effects apoe x age with ADNI interaction effects with age and apoe#
########################################################################################
dir_1 = base_directory + "/alfa_interaction_test_novol_yed/"
dir_2 = base_directory + "/adni_interaction_ageapoe_novol_yed/"
out_dir = output_directory + "/3_ALFA_int_ADNI_int/"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

dir_name_1 = "ALFA-int"
dir_name_2 = "ADNI-int"

tests = [
    ("agesqNC_HE_HO", "agesqNC_HE_HO"),
    ("agesqHE+HO-NC", "agesqHE+HO-NC")]

# Folders of the hippocampal sides
hippo_folders = [('lefthippo_ageapoe', 'lefthippo_ageapoe'), ('righthippo_ageapoe', 'righthippo_ageapoe')]

# Run the actual tests
run_tests(out_dir, hippo_folders, dir_name_1, dir_name_2, tests)


#########################################################
#4. ALFA baseline effects with ADNI NC baseline effects.#
#########################################################
dir_1 = base_directory + "/alfa_baseline_yed_novol/"
dir_2 = base_directory + "/adni_nc_baseline_test_novol_yed/"
out_dir = output_directory + "/4_ALFA_ADNI_NC_baseline/"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

dir_name_1 = "ALFA"
dir_name_2 = "ADNI-NC"

tests = [
    ("+age", "+age"),
    ("+female_-male", "+female_-male"),
    ("NC_HE_HO", "NC_HE_HO"),
    ("HE+HO-NC", "HE+HO-NC")]

# Folders of the hippocampal sides
hippo_folders = [('lefthippo', 'lefthippo'), ('righthippo', 'righthippo')]

# Run the actual tests
run_tests(out_dir, hippo_folders, dir_name_1, dir_name_2, tests)


###############################################################
#5.1 ALFA interaction effects with ADNI NC interaction effects.#
###############################################################
dir_1 = base_directory + "/alfa_interaction_test_lin_novol_yed/"
dir_2 = base_directory + "/adni_nc_interaction_test_lin_novol/"
out_dir = output_directory + "/5_1_ALFA_int_ADNI_NC_int/"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

dir_name_1 = "ALFA-int-lin"
dir_name_2 = "ADNI-NC-int-lin"

tests = [
    ("ageNC_HE_HO", "ageNC_HE_HO"),
    ("ageHE+HO-NC", "ageHE+HO-NC")]

# Folders of the hippocampal sides
hippo_folders = [('lefthippo_ageapoe', 'lefthippo_ageapoe'), ('righthippo_ageapoe', 'righthippo_ageapoe')]

# Run the actual tests
run_tests(out_dir, hippo_folders, dir_name_1, dir_name_2, tests)

###############################################################
#5.2 ALFA interaction effects with ADNI NC interaction effects.#
###############################################################
dir_1 = base_directory + "/alfa_interaction_test_novol_yed/"
dir_2 = base_directory + "/adni_nc_interaction_test_novol/"
out_dir = output_directory + "/5_2_ALFA_int_ADNI_NC_int/"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

dir_name_1 = "ALFA-int-lin"
dir_name_2 = "ADNI-NC-int-lin"

tests = [
    ("agesqNC_HE_HO", "agesqNC_HE_HO"),
    ("agesqHE+HO-NC", "agesqHE+HO-NC")]

# Folders of the hippocampal sides
hippo_folders = [('lefthippo_ageapoe', 'lefthippo_ageapoe'), ('righthippo_ageapoe', 'righthippo_ageapoe')]

# Run the actual tests
run_tests(out_dir, hippo_folders, dir_name_1, dir_name_2, tests)

#####################################################
#6. ALFA baseline effects with ALL baseline effects.#
#####################################################
dir_1 = base_directory + "/alfa_baseline_yed_novol/"
dir_2 = base_directory + "/all_baseline_test_novol_yed/"
out_dir = output_directory + "/6_ALFA_ALL_baseline/"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

dir_name_1 = "ALFA"
dir_name_2 = "all"

tests = [
    ("+age", "+age"),
    ("+female_-male", "+female_-male"),
    ("NC_HE_HO", "NC_HE_HO"),
    ("HE+HO-NC", "HE+HO-NC"),
    ("NC_HE_HO", "AD_MCI_CN")]

# Folders of the hippocampal sides
hippo_folders = [('lefthippo', 'lefthippo'),('righthippo', 'righthippo')]

# Run the actual tests
run_tests(out_dir, hippo_folders, dir_name_1, dir_name_2, tests)

###########################################################
#7. ALFA interaction effects with ALL interaction effects.#
###########################################################
dir_1 = base_directory + "/alfa_interaction_test_novol_yed/"
dir_2 = base_directory + "/all_interaction_test_ageapoe_novol/"
out_dir = output_directory + "/7_ALFA_int_ALL_int/"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

dir_name_1 = "ALFA-int"
dir_name_2 = "ALL-int"

tests = [
    ("agesqNC_HE_HO", "agesqNC_HE_HO"),
    ("agesqHE+HO-NC", "agesqHE+HO-NC")]

# Folders of the hippocampal sides
hippo_folders = [('lefthippo_ageapoe', 'lefthippo_ageapoe'), ('righthippo_ageapoe', 'righthippo_ageapoe')]

# Run the actual tests
run_tests(out_dir, hippo_folders, dir_name_1, dir_name_2, tests)


#######################################################################
#8. ALFA interaction effects with ADNI interaction effects of edat DX.#
#######################################################################

dir_1 = base_directory + "/alfa_interaction_test_novol_yed/"
dir_2 = base_directory + "/adni_interaction_test_dxage_novol_yed/"
out_dir = output_directory + "/8_ALFA_int_ADNI_dxage/"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

dir_name_1 = "ALFA-int"
dir_name_2 = "ALL-int"


tests = [
    ("agesqNC_HE_HO", "agesqCN_MCI_AD"),
    ("agesqNC_HE_HO", "agesq-CN_MCI_AD"),
    ("agesqHE+HO-NC", "agesqMCI+AD-CN")]

# Folders of the hippocampal sides
hippo_folders = [('lefthippo_ageapoe', 'lefthippo_ageapoe'), ('righthippo_ageapoe', 'righthippo_ageapoe')]

# Run the actual tests
run_tests(out_dir, hippo_folders, dir_name_1, dir_name_2, tests)
"""