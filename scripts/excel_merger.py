"""This script takes the outputs of best_model.py and generates an excel
sheet by merging the bluer and redder regimes of the same mask."""

import numpy as np
import pandas as pd 
import sys

#user input of blue and red masknames
maskname_b = sys.argv[1]
maskname_r = sys.argv[2]

#relative path to find the best_model.py outputs
relative_path = "../results/Max_z_n_width/" 
df_b = pd.read_csv(relative_path + maskname_b + ".txt")
df_r = pd.read_csv(relative_path + maskname_r + ".txt")
#
# #rename columns
df_b.columns = ["Slit Num", "RA_b", "Dec_b", "z_b", "vel_b", "amp_b"]
df_r.columns = ["Slit Num", "RA_r", "Dec_r", "z_r", "vel_r", "amp_r"]

#merge dataframes on slit number
df_combined = df_b.merge(df_r, on = "Slit Num")

#convert index to slit number
df_combined = df_combined.set_index("Slit Num")

#insert blank columns for confidence and notes
df_combined.insert(loc = 5, column = "Confidence_b", value = "")
df_combined.insert(loc = 6, column = "Notes_b", value = "")
df_combined['Confidence_r'] = ""
df_combined['Notes_r'] = ""

#replace -999 values with nothing
replace_999 = {-999:""}
df_combined = df_combined.replace(replace_999)

#drop last row which is an artifact of coding output from best_model.py
df_combined = df_combined[:-1]

#generate output file
fname = relative_path + maskname_b + "+" + maskname_r + "_visual-inspected.xlsx"
df_combined.to_excel(fname)