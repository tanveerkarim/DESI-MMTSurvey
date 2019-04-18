""""This script takes outputdata/ files as input and generates 
best models (redz and width) for a given mask"""

import numpy as np
from utils_spec1d import datareader, wave_grid, lambda_to_z
from global_var import *
import sys
import matplotlib.pyplot as plt
plt.style.use('ggplot')

def bestModel(maskname, idx, z, SNRdata, delChi2data):
	"""Returns txt file of best model described by redz and width
	per slit for a given mask
	----------
	maskname: name of the mask + '-' + grating number
	idx: index of a slit for a given maskname
	z: 0th output of the SNR_calculator function; redshift range
	SNRdata: 2nd output of the SNR_calculator function; SNR data cube
	delChi2data: 3rd output of the SNR calculator function; delChi2 data cube
					for plotting purposes
	"""	
	
	#trim to not let artifacts mess up hypothesis testing
	trimleft = (image[idx] != 0).argmax()
	trimright = np.max(np.nonzero(image[idx]))
	
	#since image size and delChi2data not same size, we need to convert
	#to z space to do the trimming
	zleft = np.abs(z - (wg[trimleft]/lambda0 - 1)).argmin()
	zright = np.abs(z - (wg[trimright]/lambda0 - 1)).argmin() 
	
	delChi2data_trimmed = delChi2data[:, :, zleft:zright]
	
	#if(idx == 55):
		#print(image[idx].shape)
		#print(delChi2data[idx].shape)
		#print(delChi2data_trimmed[idx].shape)
		#print(np.min(delChi2data[idx]))
		#print(np.min(delChi2data_trimmed[idx]))
		#print(delChi2data[idx][:, 0])
		#print(delChi2data_trimmed[idx][:, 0])
	#Find width and z indices for highest neg delta chi2 change
	w_idx, redshift_idx = np.argwhere(delChi2data[idx] == np.min(delChi2data_trimmed[idx]))[0]
	
	if(SNRdata[idx, w_idx, redshift_idx] >= 10):
		return z[redshift_idx], sigma_v[w_idx], redshift_idx, w_idx
	else:
		return -999, -999, -999, -999

#user input of maskname; as str
maskname = sys.argv[1]

#read in data based on maskname
data = datareader(maskname)
image = data['data_ivar'][:, 0, :]
datarows = len(image)

#Call following values from outputdata/: 
#z_range, widths, SNRdata, Ampdata, delChi2data
SNRdata = np.load("../results/outputdata/" + maskname + "/" + maskname + "-" + str(relative_strength) + "-SNR.npy")
Ampdata = np.load("../results/outputdata/" + maskname + "/" + maskname + "-" + str(relative_strength) + "-Amp.npy")
delChi2data = np.load("../results/outputdata/" + maskname + "/" + maskname + "-" + str(relative_strength) + "-delChi2.npy")
z_range = np.load("../results/outputdata/" + maskname + "/" + maskname + "-" + str(relative_strength) + "-z_range.npy")
widths = np.load("../results/outputdata/" + maskname + "/" + maskname + "-" + str(relative_strength) + "-widths.npy")

#Initalise arrays to store redshift and dispersion velocity values
zmax = np.zeros(datarows)
vmax = np.zeros(datarows)
redz_idx = np.zeros(datarows)
widths_idx = np.zeros(datarows)
crval1, wg = wave_grid(data) #wavelength grid of the 1d spectra

from time import time
start = time()
for i in range(1, datarows):
	zmax[i-1], vmax[i-1], redz_idx[i-1], widths_idx[i-1] = bestModel(maskname, idx = i, z=z_range,\
	SNRdata=SNRdata, delChi2data = delChi2data)
end = time()

tot_time = end - start
print(f'Total time: {tot_time}')

#Save redshift and width values in a txt file
import pandas as pd

df = pd.DataFrame({'max_z':zmax, 'max_w':vmax})
df.index += 1
df.to_csv('../results/Max_z_n_width/' + maskname + ".txt")

#save redshift and width idx for peak-zoom.py and hyp-test.py 
redz_idx = redz_idx.astype(int)
widths_idx = widths_idx.astype(int)
tmp = np.hstack((redz_idx[:, np.newaxis], widths_idx[:, np.newaxis]))
np.save("../results/Max_z_n_width/" + maskname + '-indices.npy', tmp)

#Generate initial histograms of redz and widths
plt.hist(zmax[zmax >0], bins = 15, facecolor = 'red', alpha = 0.75)
plt.xlabel('Redshift')
plt.ylabel('Frequency')
plt.title('Redshift histogram of ' + maskname)
plt.savefig('../results/histograms/redshift/' + maskname + '.png', dpi = 200, bbox_inches = None)
plt.close()

plt.hist(vmax[zmax > 0], bins = 10, facecolor = 'red', alpha = 0.75)
plt.xlabel('Width')
plt.ylabel('Frequency')
plt.title('Width histogram of ' + maskname)
plt.savefig('../results/histograms/width/' + maskname + '.png', dpi = 200, bbox_inches = None)
plt.close()