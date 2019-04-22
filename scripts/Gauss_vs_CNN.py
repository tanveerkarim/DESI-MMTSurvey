"""This script compares redshifts measured by my method (Gaussian filter) 
to the redshifts measured by Jae's method (CNN)."""

import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
plt.style.use("fivethirtyeight")

#read Jae's data; this is saved in results/JaeResults/union-catalog-results/
df_CNN = pd.read_csv("../data/JaeResults/union-catalog.csv")

#user inputs actual name of mask
maskname = sys.argv[1]
MMT_maskID_b = sys.argv[2]

#filter out relevant mask and turn it into a pandas dataframe
df_CNN = df_CNN[df_CNN['mask_name'] == maskname]
df_CNN = df_CNN[1:-1] #clip first and last point to match Gaussian_Filter data length

#turn -999 to np.nan
replace_999 = {-999:np.nan}
df_CNN = df_CNN.replace(replace_999)

#read in Gaussian filter results
df_Gauss = pd.read_excel("../results/Max_z_n_width/" + MMT_maskID_b + "+" +\
 str(int(MMT_maskID_b) + 1) + "_visual-inspected.xlsx")

#take average between blue and red mask and generate final redz
df_Gauss['z_final'] = df_Gauss[['z_b', 'z_r']].mean(axis = 1)

#redshift arrays
jz = np.array(df_CNN['redz'])
tz = np.array(df_Gauss['z_final'])

#plot comparison
x = np.arange(0, 1.7, .1)
alpha = 0.7
plt.plot(jz, tz, marker = "o", ls = "", alpha = alpha)
plt.plot(x, x, alpha = alpha)
plt.xlabel("CNN redshifts")
plt.ylabel("Gaussian redshifts")
fname = "../results/summary_statistics/CNN_vs_Gauss/" + maskname \
+ "-CNN_vs_Gauss_comparison.png"
plt.savefig(fname, dpi = 250, bbox_inches = 'tight')

#check if there are any discrepancies over 0.01 in redshift
idx = np.where(jz-tz > 0.01)
discrepant_slits_G = df_Gauss.iloc[idx]
discrepant_slits_C = df_CNN.iloc[idx]
discrepant_slits = pd.concat([discrepant_slits_G, discrepant_slits_C], axis = 1)

#show user list of discrepancies
print('discrepant slits:')
print('-----------------')
print(discrepant_slits)
print("------------------------------------------------------")

#Show stats on how many redshifts CNN method did not identify
idx = np.where(np.isnan(jz))
print("slits where Gaussian filter found redz but CNN did not:")
print("Total number of slits: ", np.sum((np.isnan(tz[np.isnan(jz)]) == False)))
print("------------------------------------------------------")
print(df_Gauss.iloc[idx].dropna(subset = ['z_final']))
print("------------------------------------------------------")

#Show stats on how many redshifts Gaussian Filter method did not identify
idx = np.where(np.isnan(tz))
print("slits where CNN found redz but Gaussian Filter did not:")
print("Total number of slits: ", np.sum((np.isnan(jz[np.isnan(tz)]) == False)))
print("------------------------------------------------------")
print(df_CNN.iloc[idx].dropna())
print("------------------------------------------------------")