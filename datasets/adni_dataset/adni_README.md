## Readme for ADNI dataset
Explanation of the different files in this directory. For all the scripts to work, paths to the required files at the start of this script should be updated with the paths of the required files.

**Exploration of the .csv data:** We used the following scripts. They do similar things, and are used mainly to gather information of the cohort.

* *adni_data_exploration.py*
* *adni_data_study.py*

**Processing of the meshes:** These scripts process and extract the meshes, separate to left and right, apply Procustes analysis and saves the resulting meshes. Dependind on the subset of the cohort used, we use different scripts: 

* *adni_preprocessing_onlyCN.py*: use only the cognitively healthy subjects.
* *adni_preprocessing_young.py*: Use only young subjects.
* *adni_preprocessing.py*: Use all subjects of the cohort.
* *adni_apoe_risk.py*: Add a column to the data with the APOE OR risk.