"""
Data study full info

We now have more information proportioned by Juando on our dataset and is probably necessary to do an extra exploration of our data.
Which features to keep? Which to remove?

"""
import pandas as pd
import os
import numpy as np
from scipy import stats
import seaborn
from statsmodels.formula.api import ols

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

def show_summary(df):
    """Show a summary of gender and age for the input dataframe."""
    df = df[["gender", "age", "education_years", "vol_r", "vol_l"]]
    df["vol_total"] = df["vol_r"] + df["vol_l"]
    print("Number of subjects: " + str(len(df)))
    print("Age: " + str(round(df.age.mean(), 2)) + "+" + str(round(df.age.std(),2)))
    print("Education years: " + str(round(df.education_years.mean(), 2)) + "+" + str(round(df.education_years.std(),2)))
    print("Females " + str(len(df[df.gender ==' F'])))
    print("Gender % (females) " + str(round(len(df[df.gender ==' F'])/len(df)*100, 2)))
    print("Bilateral volume " + str(round(df.vol_total.mean(), 2)) + "+" + str(round(df.vol_total.std(),2)))

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
    df_noncarr = df_apoe[df_apoe["APOE_e4_num"] == 0]
    show_summary(df_noncarr)

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
full_summary(df_apoe)

# Do inferential statistics on each of the covariates

# Age
model = ols('APOE_e4_num ~ age', df_apoe).fit()
print(model.summary())

# Education
model = ols('APOE_e4_num ~ education_years', df_apoe).fit()
print(model.summary())

# Volume
df_apoe["vol_total"] = df_apoe["vol_r"] + df_apoe["vol_l"]
model = ols('APOE_e4_num ~ vol_total', df_apoe).fit()
print(model.summary())

# Gender (need chi square statistic)
female_nc = len(df_apoe[(df_apoe.gender == ' F') & (df_apoe.APOE_e4_num == 0)])
female_ho = len(df_apoe[(df_apoe.gender == ' F') & (df_apoe.APOE_e4_num == 1)])
female_he = len(df_apoe[(df_apoe.gender == ' F') & (df_apoe.APOE_e4_num == 2)])

male_nc = len(df_apoe[(df_apoe.gender == ' M') & (df_apoe.APOE_e4_num == 0)])
male_ho = len(df_apoe[(df_apoe.gender == ' M') & (df_apoe.APOE_e4_num == 1)])
male_he = len(df_apoe[(df_apoe.gender == ' M') & (df_apoe.APOE_e4_num == 2)])

table = np.array([[male_nc, male_ho, male_he], [female_nc, female_ho, female_he]])

chisq, p, dof, exp = stats.chi2_contingency(table)
print(chisq)
print(p)
