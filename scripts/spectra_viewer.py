"""This script generates postage stamp given a mask, slit idx and wavelength value"""
import numpy as np
from global_var import *
from utils_spec1d import datareader, npy_reader, Window, wave_grid, Model, lambda_to_z, extract_single_data, preprocess_bino
import sys
import matplotlib.pyplot as plt
plt.style.use("fivethirtyeight")

#user input of maskname; as str
maskname = sys.argv[1]
idx = np.int(sys.argv[2])
central_wave = np.float(sys.argv[3])
input_data_dir = "../../../../DATA_MAY18/" #MAKE SURE TO CHANGE THIS DEPENDING ON WHERE 2D FILES ARE

data_err, list_headers = preprocess_bino(data_dir = input_data_dir + maskname + "/", masknumber=maskname, pos_width = pos_width)
dat, err, header = extract_single_data(data_err, list_headers, idx)

print(dat.shape)
#generate wavegrid for the given mask
data = datareader(maskname)
_, wg = wave_grid(data)

#32 pixel wide postage stamp
rang = 60

#find closest pixel to central_wave 
arg = np.abs(wg - central_wave).argmin()

#slice data according to the central_wave and range
dat_windowed = dat[:,arg-rang:arg+rang+1]
wg_windowed = wg[arg-rang:arg+rang+1]

#show plot
plt.imshow(dat_windowed, cmap = 'gray',\
						vmin = -1, vmax = 1, aspect = aspect, \
						extent = (wg_windowed[0], wg_windowed[-1], dat.shape[0], 0))
plt.grid(False)
plt.show()