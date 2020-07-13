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

alfa = ""

df_alfa = pd.read_csv(alfa)

print(df_alfa.head())

# Count the APOE_STATUS info
def apoe_summary(df_apoe, apoe_status):
    # e2e2
    print(apoe_status)
    df_short = df_apoe[df_apoe["APOE_status"] == apoe_status]
    print(len(df_short))

apoe_summary(df_alfa, "Apo-ε2/ε2")
apoe_summary(df_alfa, "Apo-ε2/ε3")
apoe_summary(df_alfa, "Apo-ε3/ε3")
apoe_summary(df_alfa, "Apo-ε2/ε4")
apoe_summary(df_alfa, "Apo-ε3/ε4")
apoe_summary(df_alfa, "Apo-ε4/ε4")

print("Total subjects:")
print(len(df_alfa))
# Create dictionary to change string to simplified allele code
dict_allele = {
    "Apo-ε2/ε2": "22",
    "Apo-ε2/ε3": "23",
    "Apo-ε3/ε3": "33",
    "Apo-ε2/ε4": "24",
    "Apo-ε3/ε4": "34",
    "Apo-ε4/ε4": "44"
}

# Add the two new columns
df_alfa['APOE_allele'] = df_alfa['APOE_status'].map(dict_allele)

# Create dictionary to convert the specific combination to the OR value
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

df_alfa['APOE_OR_old'] = df_alfa['APOE_allele'].map(dict_or)
df_alfa['APOE_OR'] = np.log(df_alfa['APOE_allele'].map(dict_or_def))

print(len(df_alfa))
