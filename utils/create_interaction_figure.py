# %% Create interaction figure
# Create interaction graphics between age (or any other covariate) and effect. This function is supposed to be called from matlab.
# WARNING: DIRTY
import numpy as np
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt

csv_mesh = sys.argv[1]
out_path = sys.argv[2]
#To recenter data again
# mean_age = float(sys.argv[3])
mean_age = 57.6639

values = pd.read_csv(csv_mesh, header=None)
values.columns = ["Y (Adj)", "Age", "Apoe", "LinMod"]

values["Y (Adj)"] = abs(values["Y (Adj)"])
values["Age"] = values["Age"]# + mean_age

#print(values["Y"])

# Així que he de guardar-ne 3 i carregar-les.

params = values.LinMod.values
print(params[0:7])
# print(params[0])
# Create the line
# Line must be created w.r.t. the non-collinerar
# x = list(range(-10, 15, 1))
x = list(range(45, 75, 1))

y_nc = [float(t)*params[1] + float(t)*float(t)*params[2] + params[3]*0.0 + float(t)*0.0*params[4] + float(t)*float(t)*0.0*params[5] for t in x]
y_he = [float(t)*params[1] + float(t)*float(t)*params[2] + params[3]*1.0 + float(t)*1.0*params[4] + float(t)*float(t)*1.0*params[5] for t in x]
y_ho = [float(t)*params[1] + float(t)*float(t)*params[2] + params[3]*2.0 + float(t)*2.0*params[4] + float(t)*float(t)*2.0*params[5] for t in x]

# First verison of the figure
sns.set(style="darkgrid")

g = sns.lmplot(x="Age", y="Y (Adj)", hue="Apoe", order=2, data=values, scatter_kws={'s':18}, fit_reg=True, truncate=True, ci=0)
ax = plt.gca()
# Only for HE NC HO
# ax.set_ylim(7.25, 10.5)
# ax.set_title("P-value: " + str(params[6]))
g.savefig(out_path, dpi=300)

plt.figure()
plt.plot(x,y_nc, c='b')
plt.plot(x,y_he, c='r')
plt.plot(x,y_ho, c='g')
plt.savefig(out_path + '_param.png', dpi=300)
