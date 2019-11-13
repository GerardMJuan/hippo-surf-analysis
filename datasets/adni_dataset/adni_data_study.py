"""
Script to study ADNI data. This is meant to be run over the .csv obtained in adni_data_exploration
"""
import pandas as pd
import os
import numpy as np
from scipy import stats
import seaborn
from statsmodels.formula.api import ols

# Add the path to the .csv file that cointains the adni csv.
ADNI_csv = ""

df_adni = pd.read_csv(ADNI_csv)
print(df_adni.columns)

# Volume
df_adni["vol_total"] = df_adni["vol_r"] + df_adni["vol_l"]


# Assumim que la selecció feta anteriorment és la correcta: és a dir,
# que hem seleccionat els meshes que volem i ara l'únic que fem és
 # agafar les columnes extra desitjades.

# Columnes a afegir (de moment aquestes), i iD, per poder fer el merge
# add_columns =  ["ID", "education_years", "APOE_status", "APOE_e4_num"]

# Fem un merge dels dos datasets a partir de la columna ID i aquests dos datasets.
# df_apoe_fused = df_apoe.merge(df_alfa_full[add_columns], how='left', on='ID')
# df_apoe_fused.head()
# Comprovem quants NA hi ha a cada una de les columnes noves
# print(len(df_apoe_fused[df_apoe_fused['education_years'] == np.nan]))

def show_summary(df):
    """Show a summary of gender and age for the input dataframe."""
    df = df[["PTGENDER", "AGE", "PTEDUCAT", "vol_r", "vol_l", "MMSE"]]
    df["vol_total"] = df["vol_r"] + df["vol_l"]
    print("Number of subjects: " + str(len(df)))
    print("Age: " + str(round(df.AGE.mean(),2)) + "+" + str(round(df.AGE.std(),2)))
    print("Education years: " + str(round(df.PTEDUCAT.mean(),2)) + "+" + str(round(df.PTEDUCAT.std(),2)))
    print("Females " + str(len(df[df.PTGENDER =='Female'])))
    print("Gender % (females) " + str(round(len(df[df.PTGENDER =='Female'])/len(df)*100, 2)))
    print("Bilateral volume " + str(round(df.vol_total.mean(),2)) + "+" + str(round(df.vol_total.std(),2)))
    print("MMSE " + str(round(df.MMSE.mean(),2)) + "+" + str(round(df.MMSE.std(),2)))

## Non carriers
def full_summary(df_apoe):
    # print("Full summary of the data")
    # # e2e2
    # print('Non-carrier, e2e2')
    # df_e2e2 = df_apoe[df_apoe["APOE_status"] == "Apo-ε2/ε2"]
    # show_summary(df_e2e2)
    #
    # # e2e3
    # print('Non-carrier, e2e3')
    # df_e2e3 = df_apoe[df_apoe["APOE_status"] == "Apo-ε2/ε3"]
    # show_summary(df_e2e3)
    #
    # # e3e3
    # print('Non-carrier, e3e3')
    # df_e3e3 = df_apoe[df_apoe["APOE_status"] == "Apo-ε3/ε3"]
    # show_summary(df_e3e3)

    # Total non carriers
    print('Total non carriers')
    df_noncarr = df_apoe[df_apoe["APOE4"] == 0]
    show_summary(df_noncarr)
    print(' ')

    ## Carriers
    # # e2e4
    # print('Heterozygotes, e2e4')
    # df_e2e4 = df_apoe[df_apoe["APOE_status"] == "Apo-ε2/ε4"]
    # show_summary(df_e2e4)
    #
    # # e3e4
    # print('Heterozygotes, e3e4')
    # df_e3e4 = df_apoe[df_apoe["APOE_status"] == "Apo-ε3/ε4"]
    # show_summary(df_e3e4)

    # Total non carriers
    print('Total Heterozygotes')
    df_noncarr = df_apoe[df_apoe["APOE4"] == 1]
    show_summary(df_noncarr)
    print(' ')

    # e4e4
    print('Homozygotes, e4e4')
    df_e4e4 = df_apoe[df_apoe["APOE4"] == 2]
    show_summary(df_e4e4)
    print(' ')

    # Total non carriers
    print('Total number of subjects')
    show_summary(df_apoe)
    print(' ')


"""
We need to do a full summary, and a summary for:
CN
MCI
AD

So that all can be included in the table in the paper.
"""
# Tornem a fer les estadistiques del article anterior.
print('###Summary of TOTAL###')
full_summary(df_adni)

# Summaries of each
print('###Summary of CN###')
full_summary(df_adni[df_adni.DX_bl == 'CN'])
print('###Summary of MCI###')
full_summary(df_adni[df_adni.DX_bl == 'LMCI'])
print('###Summary of AD###')
full_summary(df_adni[df_adni.DX_bl == 'AD'])

# Do inferential statistics on each of the covariates

# Age
model = ols('apoe_int ~ AGE', df_adni).fit()
print(model.summary())

# Education
model = ols('apoe_int ~ PTEDUCAT', df_adni).fit()
print(model.summary())

# Volume
model = ols('apoe_int ~ vol_total', df_adni).fit()
print(model.summary())

# Volume
model = ols('apoe_int ~ MMSE', df_adni).fit()
print(model.summary())

# Gender (need chi square statistic)
female_nc = len(df_adni[(df_adni.PTGENDER == 'Female') & (df_adni.APOE4 == 0)])
female_ho = len(df_adni[(df_adni.PTGENDER == 'Female') & (df_adni.APOE4 == 1)])
female_he = len(df_adni[(df_adni.PTGENDER == 'Female') & (df_adni.APOE4 == 2)])

male_nc = len(df_adni[(df_adni.PTGENDER == 'Male') & (df_adni.APOE4 == 0)])
male_ho = len(df_adni[(df_adni.PTGENDER == 'Male') & (df_adni.APOE4 == 1)])
male_he = len(df_adni[(df_adni.PTGENDER == 'Male') & (df_adni.APOE4 == 2)])

table = np.array([[male_nc, male_ho, male_he], [female_nc, female_ho, female_he]])

chisq, p, dof, exp = stats.chi2_contingency(table)
print(chisq)
print(p)
