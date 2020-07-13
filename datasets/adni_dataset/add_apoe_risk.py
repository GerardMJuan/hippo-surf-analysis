"""
Quick script that adds a column to the 
APOE data associated with the risk from the paper:

Reiman, E. M., Arboleda-Velasquez, J. F., Quiroz, Y. T., Huentelman, M. J., Beach, T. G., Caselli, R. J., … Zhao, Y. (2020). 
Exceptionally low likelihood of Alzheimer’s dementia in APOE2 homozygotes from a 5,000-person neuropathological study. Nature Communications, 11(1). 
https://doi.org/10.1038/s41467-019-14279-8

For later experimentation.

Generated on Windows
"""
import os
import pandas as pd
import numpy as np

adni = ""
adni_apoe = ""

df_adni = pd.read_csv(adni)
df_adni_apoe = pd.read_csv(adni_apoe)

add_columns = ["RID", "APGEN1", "APGEN2"]

# Fem un merge dels dos datasets a partir de la columna ID i aquests dos datasets.
df_apoe_fused = df_adni.merge(df_adni_apoe[add_columns], how='left', on='RID')
df_apoe_fused.head()

#Check unique occurences of each
print(df_apoe_fused["APGEN1"].value_counts())
print(df_apoe_fused["APGEN2"].value_counts())

df_apoe_fused["APOE_allele"] = df_apoe_fused["APGEN1"].astype(str) + df_apoe_fused["APGEN2"].astype(str)

print(df_apoe_fused["APOE_allele"].value_counts())

# Dict with 23 as reference (same as on the paper)
dict_or = {
    "22": 0.34,
    "23": 1.0,
    "33": 2.6,
    "24": 6.96,
    "34": 15.92,
    "44": 81.05
}

#DICT OR WITH 22 AS REFERENCE
dict_or_def = {
   "22": 1.0,
   "23": 2.94,
   "33": 7.65,
   "24": 20.47,
   "34": 46.82,
   "44": 238.38
}


df_apoe_fused['APOE_OR_old'] = df_apoe_fused['APOE_allele'].map(dict_or)

df_apoe_fused['APOE_OR'] = np.log(df_apoe_fused['APOE_allele'].map(dict_or_def))

print(len(df_apoe_fused))
