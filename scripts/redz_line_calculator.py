"""This script calculates redshift assuming a different line hypothesis is correct.
best-model.py always assumes that the best model is a [OII] doublet. But, based 
on hyp_test.py, we sometimes see that in fact the best-model is a different line, 
e.g. [O III] 49 or [O III] 50."""

import numpy as np
import pandas as pd
from global_var import *
from utils import datareader, npy_reader, Window, wave_grid, Model, lambda_to_z, extract_single_data, preprocess_bino
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.style.use("ggplot")

#str to values dict
species = dict({'hb':hB0, 'o1':OIII1,
						'o2':OIII2, 'ha':ha})

#user input of maskname; as str
#User inputs in sequential order: mask, blue/red flag, slit_idx, ion i.e. the line at which we want to calculate the redz, and confidence, i.e. out of 3, how sure are we about redz measurement
maskname = sys.argv[1]
flag = (sys.argv[2].lower() == 'true') #user entered input whether bluer or redder mask
confidence = sys.argv[3]
ion = sys.argv[4]
idx = np.int(sys.argv[5])

if(flag == True):
	maskname_b = maskname
	maskname_r = str(int(maskname) + 1)
else:
	maskname_r = maskname
	maskname_b = str(int(maskname) - 1)
	
#call calculated values from SNR-Amp.py for maskname slitidx
zdata, widthsdata, SNRdata, Ampdata, delChi2data = npy_reader(maskname)

#best models indices from best-model.py
bestmodels = np.load("../results/Max_z_n_width/" + maskname + "-indices.npy")
z_bestidx, w_bestidx = bestmodels[idx-1]

#generate the new redshift hypothesis based on the line
lambda_observed = lambda0 * (1 + zdata[z_bestidx])
z_hypothesis = lambda_observed / species[ion] - 1
print(np.round(z_hypothesis,6))

#relative path to find the excel_merger.py outputs
relative_path = "../results/Max_z_n_width/"
fname = relative_path + maskname_b + "+" + maskname_r + "_visual-inspected.xlsx"
df_combined = pd.read_excel("../results/Max_z_n_width/" + maskname_b + "+" + maskname_r + "_visual-inspected.xlsx")
df_combined = df_combined.set_index("Slit Num")

#replace -999 values with nothing
replace_999 = {np.nan:""}
df_combined = df_combined.replace(replace_999)

#update value in xslx file
if(flag == 'true'):
	df_combined['z_b'][idx] = z_hypothesis
	df_combined['Confidence_b'][idx] = int(confidence)
	df_combined['Notes_b'][idx] = 'w*' + ion
else:
	df_combined['z_r'][idx] = z_hypothesis
	df_combined['Confidence_r'][idx] = int(confidence)
	df_combined['Notes_r'][idx] = 'w*' + ion
	
#generate output file
df_combined.to_excel(fname)
