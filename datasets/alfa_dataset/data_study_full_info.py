"""
Data study full info

We now have more information proportioned by Juando on our dataset and is probably necessary to do an extra exploration of our data.
Which features to keep? Which to remove?

"""
import pandas as pd
import os
import numpy as np

# .csv of ALFA with volumes of each segmented hippocampus
APOE_vols_csv = ""

# xls with full ALFA info
apoe_alfa_full = ""

# out csv
out_csv = ""

df_apoe = pd.read_csv(APOE_vols_csv)
df_alfa_full = pd.read_excel(apoe_alfa_full)
df_apoe.head()

## Assumim que la selecció feta anteriorment és la correcta: és a dir, que hem seleccionat els meshes que volem i ara l'únic que fem és agafar les columnes extra desitjades.

# Hem eliminat WMhypointensities
# Columnes a afegir (de moment aquestes), i iD, per poder fer el merge
add_columns =  ["ID", "education_years", "APOE_status", "APOE_e4_num"]

# Fem un merge dels dos datasets a partir de la columna ID i aquests dos datasets.
df_apoe_fused = df_apoe.merge(df_alfa_full[add_columns], how='left', on='ID')
df_apoe_fused.head()
# Comprovem quants NA hi ha a cada una de les columnes noves
print(len(df_apoe_fused[df_apoe_fused['education_years'] == np.nan]))
# print(len(df_apoe_fused[df_apoe_fused['WMhypointensities'] == np.nan]))

df_apoe_fused = df_apoe_fused.dropna()
print(len(df_apoe_fused))

def show_summary(df):
    """Show a summary of gender and age for the input dataframe."""
    df = df[["gender", "age"]]
    print("Number of subjects: " + str(len(df)))
    print("Age: " + str(round(df.age.mean(),1)) + "+" + str(round(df.age.std(),2)))
    print("Gender % (females) " + str(round(len(df[df.gender==' F'])/len(df)*100, 2)))

## Non carriers
def full_summary(df_apoe):
    print("Full summary of the data")
    # e2e2
    print('Non-carrier, e2e2')
    df_e2e2 = df_apoe[df_apoe["APOE_status"] == "Apo-ε2/ε2"]
    show_summary(df_e2e2)

    # e2e3
    print('Non-carrier, e2e3')
    df_e2e3 = df_apoe[df_apoe["APOE_status"] == "Apo-ε2/ε3"]
    show_summary(df_e2e3)

    # e3e3
    print('Non-carrier, e3e3')
    df_e3e3 = df_apoe[df_apoe["APOE_status"] == "Apo-ε3/ε3"]
    show_summary(df_e3e3)

    # Total non carriers
    print('Total non carriers')
    df_noncarr = df_apoe[df_apoe["APOE_e4_num"] == 0]
    show_summary(df_noncarr)

    ## Carriers
    # e2e4
    print('Heterozygotes, e2e4')
    df_e2e4 = df_apoe[df_apoe["APOE_status"] == "Apo-ε2/ε4"]
    show_summary(df_e2e4)

    # e3e4
    print('Heterozygotes, e3e4')
    df_e3e4 = df_apoe[df_apoe["APOE_status"] == "Apo-ε3/ε4"]
    show_summary(df_e3e4)

    # Total non carriers
    print('Total Heterozygotes')
    df_noncarr = df_apoe[df_apoe["APOE_e4_num"] == 1]
    show_summary(df_noncarr)


    # e4e4
    print('Homozygotes, e4e4')
    df_e4e4 = df_apoe[df_apoe["APOE_e4_num"] == 2]
    show_summary(df_e4e4)

    # Total non carriers
    print('Total number of subjects')
    show_summary(df_apoe)

# Tornem a fer les estadistiques del article anterior.
full_summary(df_apoe_fused)
df_apoe_fused.head()
# Need to change APOE_int and APOE_cat
# Just as a patch, when reprocessing all the meshes do as
n_apoe = []
c_apoe = []
for x in df_apoe_fused.apoe.values:
    print(x)
    if x in [22, 23, 24, 33]:
        n_apoe.append(0)
        c_apoe.append('NC')
    elif x == 34:
        n_apoe.append(1)
        c_apoe.append('HE')
    elif x == 44:
        n_apoe.append(2)
        c_apoe.append('HO')
    else:
        print('error?')

df_apoe_fused["apoe_int"] = n_apoe
df_apoe_fused["apoe_cat"] = c_apoe

df_apoe_fused.to_csv(out_csv, index=False, index_label=False)
