"""This script compares redshifts measured by my method (Gaussian filter) 
to the redshifts measured by Jae's method (CNN)."""

import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
plt.style.use("fivethirtyeight")

#read Jae's data; this is saved in results/JaeResults/union-catalog-results/
data_jae = np.load("../results/JaeResults/union-catalog-results.npy").item()

#user inputs actual name of mask; format available in data_jae['MASK_NAMES']
maskname = sys.argv[1]
MMT_maskID_b = sys.argv[1] #refers to the four digit code of the bluemask; can be found in docs/folder_to_mask.txt

#find mask's number in Jae's data file
for i in range(len(data_jae['MASK_NAMES'])):
    if(data_jae['MASK_NAMES'][i] == maskID):
        mask_num = i #mask_num corresponds to maskID in Jae's file
        break

#filter out relevant mask and turn it into a pandas dataframe
jae_df = pd.DataFrame(data_jae['REDZ'][data_jae['MASK_NUM'] == mask_num], columns=['jae_z'])
jae_df = jae_df[1:-1] #clip first and last point to match Gaussian_Filter data length

#turn -999 to np.nan
replace_999 = {-999:np.nan}
jae_df = jae_df.replace(replace_999)

#read in Gaussian filter results
tanveer_df = pd.read_excel("../results/Max_z_n_width/" + MMT_maskID_b + "+" +\
 str(int(MMT_maskID_b) + 1) + "_visual-inspected.xlsx")

#take average between blue and red mask and generate final redz
tanveer_df['z_final'] = tanveer_df[['z_b', 'z_r']].mean(axis = 1)

#redshift arrays
jz = np.array(jae_df['jae_z'])
tz = np.array(tanveer_df['z_final'])

#plot comparison
x = np.arange(0, 1.7, .1)
alpha = 0.7
plt.plot(jz, tz, marker = "o", ls = "", alpha = alpha)
plt.plot(x, x, alpha = alpha)
plt.xlabel("CNN redshifts")
plt.ylabel("Gaussian redshifts")
fname = "../results/summary_statistics/" + maskname + "CNN_vs_Gauss_comparison.png"
plt.plot(fname, dpi = 250, bbox_inches = 'tight')

#check if there are any discrepancies over 0.01 in redshift
idx = np.where(jz-tz > 0.01)
discrepant_slits = tanveer_df.iloc[idx]

#show user list of discrepancies
print(discrepant_slits)