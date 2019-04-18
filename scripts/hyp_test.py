"""This script generates all the zoomed in hypothesis test plots
based on the best-model.py output. Currently, it tests the following
hypothesis: [O II], Hbeta, OIII 48, OII 50 and Halpha."""

import numpy as np
from global_var import *
from utils_spec1d import datareader, npy_reader, Window, wave_grid, Model, lambda_to_z, extract_single_data, preprocess_bino
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.style.use("ggplot")

def hypothesis_tester(input_data_dir, maskname, slit_idx, z, widths, \
				image, ivar, wg, redshift_idx, w_idx, rel_str, fudge_fac, maskred_flag):
		"""Returns zoomed-in 1d spectra and inverse variance 
		plots around the maxima redshift.
		Args:
			maskname: name of the mask
			slit_idx: idx of slit in mask maskname
			z: 0th output of SNR_calculator(); redshift range
			widths: 1st output of SNR_calculator(); width range
			SNRdata: 2nd output of SNR_calculator(); SNR matrix
			Ampdata: 3rd output of SNR_calculator(); Amplitude matrix
			medians: 4th output of SNR_calculator(); median matrix
			image: 1d spectra for slit_idx
			ivar: 1d ivar for slit_idx
			wg: entire wavelength grid of image
			redshift_idx: redshift idx of best models
			w_idx: width idx of best models
			maskred_flag: whether a red mask of maskname exists; if so, check that
							for plotting purposes
		"""
		
		wg_windowed = Window(z[redshift_idx], wg)
		
		minidx = np.where(wg == np.min(wg_windowed))[0][0] 
		maxidx = np.where(wg == np.max(wg_windowed))[0][0]
			
		#Windowed wavelength grid over rang_idx
		image_windowed = np.copy(image[slit_idx][minidx:maxidx+1])
		ivar_windowed = np.copy(ivar[slit_idx][minidx:maxidx+1])
		
		######
		#Code to produce 2D spectra 
		if(maskname == '1626'):
			data_err, list_headers = preprocess_bino(masknumber=maskname, pos_width = pos_width_1626)
		else:
			data_err, list_headers = preprocess_bino(data_dir = input_data_dir + maskname + "/", masknumber=maskname, pos_width = pos_width)
		dat, err, header = extract_single_data(data_err, list_headers, idx)
		if(maskred_flag == True):
			redmask = np.str(np.int(maskname) + 1) #red mask is one up from blue mask
			data_err, list_headers = preprocess_bino(data_dir = input_data_dir + str(int(maskname) + 1) + "/", masknumber=redmask, pos_width = pos_width)
			dat_red, _, _ = extract_single_data(data_err, list_headers, idx)
			
			tmpdata = datareader(redmask)
			crval_red, wg_red = wave_grid(tmpdata) 
			cdelt_red = tmpdata['headers'][1]['CDELT1']*10 #convert to angstrom
		else: #if redder mask is default
			bluemask = np.str(np.int(maskname) - 1) #blue mask is one down from red mask

			data_err, list_headers = preprocess_bino(data_dir = input_data_dir + str(int(maskname) -1) + "/", masknumber=bluemask, pos_width = pos_width)
			dat_blue, _, _ = extract_single_data(data_err, list_headers, idx)
			
			tmpdata = datareader(bluemask)
			crval_blue, wg_blue = wave_grid(tmpdata) 
			cdelt_blue = tmpdata['headers'][1]['CDELT1']*10 #convert to angstrom
		
		#dat_windowed = dat[:, rang_idx]
		dat_windowed = np.copy(dat[:, minidx:maxidx+1])
		dat_windowed = dat_windowed.reshape((dat.shape[0],wg_windowed.shape[0]))
		######
		
		###--------------###
		#PLOTTING PART OF THE CODE
		
		maskname = str(maskname)
		
		##Check for hbeta, and OIII lines in 2D spectra##
		rang = np.median(np.arange(0, wg_windowed.shape[0]))
		rang = np.int(rang) # rang is in px space
		
		#Plot figure with subplots of different sizes
		tmpdata = datareader(maskname)
		cdelt = tmpdata['headers'][1]['CDELT1']*10 #convert to angstrom; 
				#cdelt*rang because rang is in px 
				#space so convert px to wavelength space and crdel1 is the 
				#conversion factor
		
		fig = plt.figure(1, figsize=(15,15))
		# set up subplot grid
		subgrid_plot_size = (5,5)
		gridspec.GridSpec(5,5)
		aspect = 'auto'
		
		def plot_hypothesis(emission_line, z_hypothesis, grid_loc, maskred_flag = maskred_flag):
			species = dict({hB0:r'H$\beta$', OIII1:'[OIII] 49',
						OIII2:'[OIII] 50', ha:r'H$\alpha$'})
				
			observed_wavelength = emission_line*(1 + z_hypothesis)
			
			if(maskred_flag): #if testing for bluer mask
				#Case 1: Check if data exists in the blue mask exclusively
				if(((emission_line + cdelt*rang)*(1 + z_hypothesis) > wg[0]) & \
				((emission_line + cdelt*rang)*(1 + z_hypothesis) < wg_red[0])): #make sure rightmost point exists within the wg
					idx_observed = np.abs(wg - observed_wavelength).argmin()
					dat_emission_line = dat[:, idx_observed - rang: idx_observed + rang + 1]
						
					if(np.abs(np.sum(dat_emission_line)) > 0): #if actual data exists; abs because data maybe negative
						plt.subplot2grid(subgrid_plot_size, grid_loc)
						plt.imshow(dat_emission_line, cmap = 'gray',aspect = aspect,\
								vmin = vmin, vmax = vmax, \
								extent = (wg[idx_observed - rang], \
										wg[idx_observed + rang], dat.shape[0], 0))
						plt.grid(False)
						
						if(grid_loc[0] == 0):
							plt.title(species[emission_line], fontsize = fontsize)
						
					else:
						if(grid_loc[0] == 0):
							plt.subplot2grid(subgrid_plot_size, grid_loc)
							plt.axis('off')
							plt.title(species[emission_line], fontsize = fontsize)
				
				#Case 2: Check if data exists in the red mask exclusively
				elif(((emission_line + cdelt_red*rang)*(1 + z_hypothesis) > wg[-1]) & \
					((emission_line + cdelt_red*rang)*(1 + z_hypothesis) < wg_red[-1])):
					idx_observed = np.abs(wg_red - observed_wavelength).argmin()
					dat_emission_line = dat_red[:, idx_observed - rang: idx_observed + rang + 1]
				
					if(np.abs(np.sum(dat_emission_line)) > 0): #if actual data exists; abs b/c data maybe negative
						plt.subplot2grid(subgrid_plot_size, grid_loc)
						plt.imshow(dat_emission_line, cmap = 'gray',aspect = aspect,\
								vmin = vmin, vmax = vmax, \
								extent = (wg_red[idx_observed - rang], \
										wg_red[idx_observed + rang], dat_red.shape[0], 0))
						plt.grid(False)
						
						if(grid_loc[0] == 0):
							plt.title(species[emission_line], fontsize = fontsize)
							
					else:
						if(grid_loc[0] == 0):
							plt.subplot2grid(subgrid_plot_size, grid_loc)
							plt.axis('off')
							plt.title(species[emission_line], fontsize = fontsize)
				
				#Case 3: Check the overlapping region between bluemask and redmask
				#Here, we need to check if data exists in our desired range in both 
				#blue and red; if exist in both, we arbitrarily choose to show plot
				#from redmask
				elif(((emission_line + cdelt*rang)*(1 + z_hypothesis) > wg_red[0]) & \
					((emission_line + cdelt*rang)*(1 + z_hypothesis) < wg[-1])):
					
					#Test to see if data exists in redmask
					idx_observed_r = np.abs(wg_red - observed_wavelength).argmin()
					dat_emission_line_r = dat_red[:, idx_observed_r - rang: idx_observed_r + rang + 1]
					
					#data in bluemask
					idx_observed = np.abs(wg - observed_wavelength).argmin()
					dat_emission_line = dat[:, idx_observed - rang: idx_observed + rang + 1]
					
					#condition to check data in redmask
					if(np.abs(np.sum(dat_emission_line_r)) > 0): #if actual data exists; abs because data maybe negative
						plt.subplot2grid(subgrid_plot_size, grid_loc)
						plt.imshow(dat_emission_line_r, cmap = 'gray',aspect = aspect,\
								vmin = vmin, vmax = vmax, \
								extent = (wg_red[idx_observed_r - rang], \
										wg_red[idx_observed_r + rang], dat_red.shape[0], 0))
						plt.grid(False)
					
					#If data does not exist in redmask, check to see if data exists in bluemask
					elif(np.abs(np.sum(dat_emission_line)) > 0): #if actual data exists; abs because data maybe negative
						plt.subplot2grid(subgrid_plot_size, grid_loc)
						plt.imshow(dat_emission_line, cmap = 'gray',aspect = aspect,\
								vmin = vmin, vmax = vmax, \
								extent = (wg[idx_observed - rang], \
										wg[idx_observed + rang], dat.shape[0], 0))
						plt.grid(False)
						
					else:
						if(grid_loc[0] == 0):
							plt.subplot2grid(subgrid_plot_size, grid_loc)
							plt.axis('off')
							plt.title(species[emission_line], fontsize = fontsize)
						
				else:
					if(grid_loc[0] == 0):
							plt.subplot2grid(subgrid_plot_size, grid_loc)
							plt.axis('off')
							plt.title(species[emission_line], fontsize = fontsize)
							
			else: #for red masks
					
				#Case 1: Check if data exists in the red mask exclusively
				if(((emission_line + cdelt*rang)*(1 + z_hypothesis) > wg_blue[-1]) & \
				((emission_line + cdelt*rang)*(1 + z_hypothesis) < wg[-1])): #make sure rightmost point exists within the wg
					idx_observed = np.abs(wg - observed_wavelength).argmin()
					dat_emission_line = dat[:, idx_observed - rang: idx_observed + rang + 1]
						
					if(np.abs(np.sum(dat_emission_line)) > 0): #if actual data exists; abs because data maybe negative
						plt.subplot2grid(subgrid_plot_size, grid_loc)
						plt.imshow(dat_emission_line, cmap = 'gray',aspect = aspect,\
								vmin = vmin, vmax = vmax, \
								extent = (wg[idx_observed - rang], \
										wg[idx_observed + rang], dat.shape[0], 0))
						plt.grid(False)
						
						if(grid_loc[0] == 0):
							plt.title(species[emission_line], fontsize = fontsize)
						
					else:
						if(grid_loc[0] == 0):
							plt.subplot2grid(subgrid_plot_size, grid_loc)
							plt.axis('off')
							plt.title(species[emission_line], fontsize = fontsize)
				
				#Case 2: Check if data exists in the blue mask exclusively
				elif(((emission_line + cdelt_blue*rang)*(1 + z_hypothesis) > wg_blue[0]) & \
					((emission_line + cdelt_blue*rang)*(1 + z_hypothesis) < wg[0])):
					idx_observed = np.abs(wg_blue - observed_wavelength).argmin()
					dat_emission_line = dat_blue[:, idx_observed - rang: idx_observed + rang + 1]
				
					if(np.abs(np.sum(dat_emission_line)) > 0): #if actual data exists; abs b/c data maybe negative
						plt.subplot2grid(subgrid_plot_size, grid_loc)
						plt.imshow(dat_emission_line, cmap = 'gray',aspect = aspect,\
								vmin = vmin, vmax = vmax, \
								extent = (wg_blue[idx_observed - rang], \
										wg_blue[idx_observed + rang], dat_blue.shape[0], 0))
						plt.grid(False)
						
						if(grid_loc[0] == 0):
							plt.title(species[emission_line], fontsize = fontsize)
							
					else:
						if(grid_loc[0] == 0):
							plt.subplot2grid(subgrid_plot_size, grid_loc)
							plt.axis('off')
							plt.title(species[emission_line], fontsize = fontsize)
				
				#Case 3: Check the overlapping region between bluemask and redmask
				#Here, we need to check if data exists in our desired range in both 
				#blue and red; if exist in both, we arbitrarily choose to show plot
				#from redmask
				elif(((emission_line + cdelt*rang)*(1 + z_hypothesis) > wg[0]) & \
					((emission_line + cdelt*rang)*(1 + z_hypothesis) < wg_blue[-1])):
						
					#Test to see if data exists in redmask
					idx_observed_r = np.abs(wg - observed_wavelength).argmin()
					dat_emission_line_r = dat[:, idx_observed_r - rang: idx_observed_r + rang + 1]	
					
					#data in bluemask
					idx_observed = np.abs(wg_blue - observed_wavelength).argmin()
					dat_emission_line = dat_blue[:, idx_observed - rang: idx_observed + rang + 1]
					
					if(np.abs(np.sum(dat_emission_line_r)) > 0): #if actual data exists; abs because data maybe negative
						plt.subplot2grid(subgrid_plot_size, grid_loc)
						plt.imshow(dat_emission_line_r, cmap = 'gray',aspect = aspect,\
								vmin = vmin, vmax = vmax, \
								extent = (wg[idx_observed_r - rang], \
										wg[idx_observed_r + rang], dat.shape[0], 0))
						plt.grid(False)
						
					#If data does not exist in redmask, check to see if data exists in bluemask
					elif(np.abs(np.sum(dat_emission_line)) > 0): #if actual data exists; abs because data maybe negative
						plt.subplot2grid(subgrid_plot_size, grid_loc)
						plt.imshow(dat_emission_line, cmap = 'gray',aspect = aspect,\
								vmin = vmin, vmax = vmax, \
								extent = (wg_blue[idx_observed - rang], \
										wg_blue[idx_observed + rang], dat.shape[0], 0))
						plt.grid(False)
						
					else:
						if(grid_loc[0] == 0):
							plt.subplot2grid(subgrid_plot_size, grid_loc)
							plt.axis('off')
							plt.title(species[emission_line], fontsize = fontsize)
						
				else:
					if(grid_loc[0] == 0):
							plt.subplot2grid(subgrid_plot_size, grid_loc)
							plt.axis('off')
							plt.title(species[emission_line], fontsize = fontsize)
		
		##--HYPOTHESIS: [OII]--##
		# OII
		plt.subplot2grid(subgrid_plot_size, (0,0))
		plt.title('[O II]', fontsize = fontsize)
		plt.imshow(dat_windowed, cmap = 'gray', aspect = aspect,\
						vmin = vmin, vmax = vmax, \
						extent = (wg_windowed[0], wg_windowed[-1], dat.shape[0], 0))
		plt.grid(False)
		
		# OIII1
		plot_hypothesis(OIII1, z[redshift_idx], (0,2))
		plt.yticks([])
		
		#OIII2
		plot_hypothesis(OIII2, z[redshift_idx], (0,3))
		plt.yticks([])
		
		# hB
		plot_hypothesis(hB0, z[redshift_idx], (0,1))
		plt.yticks([])
		
		# ha
		plot_hypothesis(ha, z[redshift_idx], (0,4))
		plt.yticks([])
			
		##--HYPOTHESIS: [Hb]--##
		#define new z hypothesis, z_Hb
		
		lambda_observed = np.median(wg_windowed)
		z_hypothesis = lambda_observed/hB0 - 1
		
		# Hb
		plt.subplot2grid(subgrid_plot_size, (1,1))
		plt.imshow(dat_windowed, cmap = 'gray', aspect = aspect, \
						vmin = vmin, vmax = vmax, \
						extent = (wg_windowed[0], wg_windowed[-1], dat.shape[0], 0))
		plt.yticks([])
		plt.grid(False)
		
		# OIII1
		plot_hypothesis(OIII1, z_hypothesis, (1,2))
		plt.yticks([])
			
		#OIII2
		plot_hypothesis(OIII2, z_hypothesis, (1,3))
		plt.yticks([])
		
		# OII
		plot_hypothesis(OII, z_hypothesis, (1,0))
		
		# ha
		plot_hypothesis(ha, z_hypothesis, (1,4))
		plt.yticks([])
		
		##--HYPOTHESIS: [OIII1]--##
		#define new z hypothesis, z_OIII1
		
		lambda_observed = np.median(wg_windowed)
		z_hypothesis = lambda_observed/OIII1 - 1
		
		# OIII1
		plt.subplot2grid(subgrid_plot_size, (2,2))
		plt.imshow(dat_windowed, cmap = 'gray', aspect = aspect, \
						vmin = vmin, vmax = vmax, \
						extent = (wg_windowed[0], wg_windowed[-1], dat.shape[0], 0))
		plt.yticks([])
		plt.grid(False)
		
		# Hbeta
		plot_hypothesis(hB0, z_hypothesis, (2,1))
		plt.yticks([])
			
		#OIII2
		plot_hypothesis(OIII2, z_hypothesis, (2,3))
		plt.yticks([])
		
		# OII
		plot_hypothesis(OII, z_hypothesis, (2,0))
		
		# ha
		plot_hypothesis(ha, z_hypothesis, (2,4))
		plt.yticks([])
		
		##--HYPOTHESIS: [OIII2]--##
		#define new z hypothesis, z_OIII2
		
		lambda_observed = np.median(wg_windowed)
		z_hypothesis = lambda_observed/OIII2 - 1
		
		# OIII2
		plt.subplot2grid(subgrid_plot_size, (3,3))
		plt.imshow(dat_windowed, cmap = 'gray', \
						vmin = vmin, vmax = vmax, aspect = aspect, \
						extent = (wg_windowed[0], wg_windowed[-1], dat.shape[0], 0))
		plt.yticks([])
		plt.grid(False)
		
		# Hbeta      
		plot_hypothesis(hB0, z_hypothesis, (3,1))
		plt.yticks([])
			
		#OIII1
		plot_hypothesis(OIII1, z_hypothesis, (3,2))
		plt.yticks([])
		
		# OII
		plot_hypothesis(OII, z_hypothesis, (3,0))
		
		# ha
		plot_hypothesis(ha, z_hypothesis, (3,4))
		plt.yticks([])
		
		##--HYPOTHESIS: h alpha--##
		#define new z hypothesis, z_OIII2
		
		lambda_observed = np.median(wg_windowed)
		z_hypothesis = lambda_observed/ha - 1
		
		# ha
		plt.subplot2grid(subgrid_plot_size, (4,4))
		plt.imshow(dat_windowed, cmap = 'gray',\
						vmin = vmin, vmax = vmax, aspect = aspect, \
						extent = (wg_windowed[0], wg_windowed[-1], dat.shape[0], 0))
		plt.yticks([])
		plt.grid(False)
		
		# Hbeta   
		plot_hypothesis(hB0, z_hypothesis, (4,1))
		plt.yticks([])
			
		#OIII1
		plot_hypothesis(OIII1, z_hypothesis, (4,2))
		plt.yticks([])
		
		# OII
		plot_hypothesis(OII, z_hypothesis, (4,0))
		
		# OIII2
		plot_hypothesis(OIII2, z_hypothesis, (4,3))
		plt.yticks([])
		
		
		plt.suptitle(f'Mask: {maskname} Slit: {slit_idx} z : {np.round(z[redshift_idx],2)} vel: {sigma_v[w_idx]} km/s', \
						fontsize = fontsize, fontname = 'serif')
		plt.tight_layout(rect=[0, 0.03, 1, 0.95])
		plt.subplots_adjust(wspace=0.1, hspace=0.1)
		fig_name = '../results/hyp_test/' + str(maskname) + '/' + str(maskname) + "-Slit-" \
		+ str(slit_idx) + '.png'
		fig.savefig(fig_name, dpi = 250, bbox_inches = 'tight')
		plt.close()
		#plt.show()

#user input of maskname; as str
maskname = sys.argv[1]
flag = (sys.argv[2].lower() == 'true') #user entered input whether bluer or redder mask
input_data_dir = "../../../../DATA_MAY18/" #MAKE SURE TO CHANGE THIS DEPENDING ON WHERE 2D FILES ARE

#read in data based on maskname
data = datareader(maskname)
image = data['data_ivar'][:, 0, :]
ivar = data['data_ivar'][:, 1, :]
datarows = len(image)

zdata, widthsdata, SNRdata, Ampdata, delChi2data = npy_reader(maskname)

#best models indices from best-model.py
bestmodels = np.load("../results/Max_z_n_width/" + maskname + "-indices.npy")

#generate wavegrid for the given mask
_, wg = wave_grid(data)

#check if directory to save plots exist; if not, create one
output_dir = '../results/hyp_test/' + str(maskname) + '/'
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

from time import time

start = time()
for idx in range(1, datarows):
#idx = 160
	redshift_idx = bestmodels[idx-1, 0]
	w_idx = bestmodels[idx-1, 1]
	
	if(redshift_idx != -999): #ignore nans
		if(SNRdata[idx, w_idx, redshift_idx] >= 10):
			hypothesis_tester(input_data_dir = input_data_dir, maskname = maskname, \
			slit_idx = idx, z = zdata, widths = widthsdata, \
			image = image, ivar = ivar, wg = wg, redshift_idx = redshift_idx, \
			w_idx = w_idx, rel_str = relative_strength, fudge_fac = fudge_factor, \
			maskred_flag = flag)
	 
end = time()
print(f'Time: {end - start}')