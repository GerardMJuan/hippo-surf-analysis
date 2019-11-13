"""
Explore distribution of subjects in the cohort

"""

# %%
import pandas as pd

# .csv indicating which subjects we have in the cohort. Subset of ADNIMERGE, with a PTID column.
Hipp_subj_path=""

# Path to ADNI adnimerge
ADNIMERGE_path=""

# Load the tables
df_real = pd.read_csv(Hipp_subj_path)
df_adni = pd.read_csv(ADNIMERGE_path)

# We only want to merge the column with the POE info, so
df_adni = df_adni[df_adni.VISCODE == 'bl']

# We want the following columns (to have the same info as in the alfa dataset)
df_adnimerge_apoe = df_adni[["PTID", "APOE4", "DX_bl", "AGE", "PTEDUCAT", "PTGENDER". "MMSE"]]

# And now merge it, assuming that PTID and Subject columns contain the same information.
df_adni_total = df_real[["PTID"]].merge(df_adnimerge_apoe, how="inner",on="PTID")


# %%
df_adni_total.drop_duplicates(inplace=True)

def show_summary(df):
    """Show a summary of gender, age and APOE for the input dataframe."""
    print("Number of subjects: " + str(len(df)))
    print("Age mean: " + str(round(df.AGE.mean(),1)) + "+" + str(round(df.AGE.std(),2)))
    print("Females: " + str(round(len(df[df['PTGENDER']=='Female'])/len(df),2)))

    df_apoe0 = df[df["APOE4"] == 0]
    df_apoe1 = df[df["APOE4"] == 1]
    df_apoe2 = df[df["APOE4"] == 2]
    print('Non-carriers: ' + str(len(df_apoe0)))
    print('Heterozygotes: ' + str(len(df_apoe1)))
    print('Homozygotes: ' + str(len(df_apoe2)))

# How many different DX do we have?

df_adni_total = df_adni_total.replace(['SMC', 'EMCI', 'LMCI'], ['CN','MCI','MCI'])

list_dx = set(df_adni_total.DX_bl.values)
print(list_dx)

for dx in list_dx:
    print(dx)
    # e2e2
    df_dx = df_adni_total[df_adni_total["DX_bl"] == dx]
    show_summary(df_dx)

n_apoe = []
c_apoe = []
for x in df_adni_total.APOE4.values:
    if x == 0:
        n_apoe.append(0)
        c_apoe.append('NC')
    elif x == 1:
        n_apoe.append(1)
        c_apoe.append('HE')
    elif x == 2:
        n_apoe.append(2)
        c_apoe.append('HO')
    else:
        print('error?')

df_adni_total["apoe_int"] = n_apoe
df_adni_total["apoe_cat"] = c_apoe

# Save df_adni to disk
out_csv = ""
df_adni_total.to_csv(out_csv, index=False, index_label=False)
