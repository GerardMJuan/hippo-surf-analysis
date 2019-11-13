## Readme for ALFA dataset
Explanation of the different files in this directory. For all the scripts to work, paths to the required files at the start of this script should be updated with the paths of the required files.

**Exploration of the .csv data:** We used the following scripts. They do similar things, and are used mainly to gather information of the cohort.

* *data_study_alfa.py* 
* *data_study_full_info.py*

**Processing of the meshes:** These scripts process and extract the meshes, separate to left and right, apply Procustes analysis and saves the resulting meshes. Depending on the subset of the cohort used, we use different scripts: 

* *mesh_processing.py*