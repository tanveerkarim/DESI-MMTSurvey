"""This script generates all the zoomed in plots of the best-model.py output"""

import numpy as np
from global_var import *
from utils_spec1d import datareader, npy_reader, Window, wave_grid, Model, lambda_to_z, extract_single_data, preprocess_bino
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.style.use("ggplot")

def PeakZoom(maskname, slit_idx, z, widths, SNRdata, Ampdata, delChi2data, \
				image, ivar, wg, redshift_idx, w_idx, rel_str, fudge_fac, maskred_flag):
		"""Returns zoomed-in 1d spectra and inverse variance 
		plots around the maxima redshift.
		Args:
			maskname: name of the mask
			slit_idx: idx of slit in mask maskname
			z: 0th output of SNR_calculator(); redshift range
			widths: 1st output of SNR_calculator(); width range
			Ampdata: 3rd output of SNR_calculator(); Amplitude matrix
			medians: 4th output of SNR_calculator(); median matrix
			image: 1d spectra for slit_idx
			ivar: 1d ivar for slit_idx
			wg: entire wavelength grid of image
			redshift_idx: redshift idx of best models
			w_idx: width idx of best models
			rel_str: relative strength of [O II] 27 wrt 29
			fudge_fac: fudge factor for sigma_slit
			maskred_flag: whether a red mask of maskname exists; if so, check that 
					for plotting purposes
		"""
		
		#range over which to be plotted
		"""rang_idx = np.where(np.logical_and(wg > (1+z[redshift_idx])*\
							(lambda0 - rang), \
							wg < (1+z[redshift_idx])*(lambda0 + rang)))"""
		
		wg_windowed = Window(z[redshift_idx], wg)
		
		minidx = np.where(wg == np.min(wg_windowed))[0][0] 
		maxidx = np.where(wg == np.max(wg_windowed))[0][0]
		delChi2_windowed = np.copy(delChi2data[slit_idx, w_idx])
			
		#Windowed wavelength grid over rang_idx
		image_windowed = np.copy(image[slit_idx][minidx:maxidx+1])
		ivar_windowed = np.copy(ivar[slit_idx][minidx:maxidx+1])
		
		######
		#Code to produce 2D spectra
		if(maskname == '1626'):
			data_err, list_headers = preprocess_bino(masknumber=maskname, pos_width = pos_width_1626)
		else:
			data_err, list_headers = preprocess_bino(masknumber=maskname)
		dat, err, header = extract_single_data(data_err, list_headers, idx)
		if(maskred_flag == True):
			redmask = np.str(np.int(maskname) + 1) #red mask is one up from blue mask
			data_err, list_headers = preprocess_bino(masknumber=redmask)
			dat_red, _, _ = extract_single_data(data_err, list_headers, idx)
			
			tmpdata = datareader(redmask)
			crval_red, wg_red = wave_grid(tmpdata) 
			cdelt_red = tmpdata['headers'][1]['CDELT1']*10 #convert to angstrom
		
		dat_windowed = np.copy(dat[:, minidx:maxidx+1])
		dat_windowed = dat_windowed.reshape((dat.shape[0],wg_windowed.shape[0]))
		######
		
		#model based on the highest SNR
		model = Model(z[redshift_idx], wg_windowed, widths[redshift_idx][w_idx], \
							Amp=Ampdata[idx, w_idx, redshift_idx], \
										relative_strength = rel_str)[:,0]
		
		#plot delchi2#
		z_grid = lambda_to_z(wg) #Convert wavelength space to redshift space to set z_range
		z_range = np.arange(z_grid[0], z_grid[-1], 0.0002) #Set redshift hypothesis grid
		
		leftidx = np.abs(wg_windowed[0] - lambda0*(1 + z_range)).argmin()
		rightidx = np.abs(wg_windowed[-1] - lambda0*(1 + z_range)).argmin()
		delChi2data = delChi2data[:, :, leftidx:rightidx+1]
		
		z_range = z_range[leftidx:rightidx+1]
		##
		
		###--------------###
		#PLOTTING PART OF THE CODE
		
		maskname = str(maskname)
		
		
		##Check for h-beta, and OIII lines in 2D spectra##
		crdelt1 = 0.62 #ad hoc value crdelt1*rang because rang is in px 
            #space so convert px to wavelength space and crdel1 is the 
            #conversion factor
    
		"""set range that will plot same window width as that of [O II] doublet.
		given the index of a given line, we will generate 2D plots of 
		xlim = [idx - crdelt1*rang, idx + crdelt1*rang]"""
		rang = np.median(np.arange(0, wg_windowed.shape[0]))
		rang = np.int(rang)
		
		#Plot figure with subplots of different sizes
		vmin = -0.5
		vmax = 1.
		
		fig = plt.figure(1)
		# set up subplot grid
		gridspec.GridSpec(4,4)
		
		# OII
		plt.subplot2grid((4,4), (0,0))
		plt.title('[O II]')
		plt.imshow(dat_windowed, cmap = 'gray', aspect = 'auto',\
						vmin = vmin, vmax = vmax, \
						extent = (wg_windowed[0], wg_windowed[-1], dat.shape[0], 0))
		plt.grid(False)
		#plt.axvline(lambda0 * (1 + z[redshift_idx]), c = 'red', linestyle = "--", alpha = alpha)
		
		# hB
		observed_wavelength = hB0*(1 + z[redshift_idx])
		
		plt.subplot2grid((4,4), (0,3))
		plt.title(r'H$\beta$')
		
		if(maskred_flag): #if bluer mask 
			if((hB0 + crdelt1*rang)*(1 + z[redshift_idx]) < wg[-1]): #make sure rightmost point exists within the wg
				idx_observed = np.abs(wg - observed_wavelength).argmin()
				dat_hB0 = dat[:, idx_observed - rang: idx_observed + rang + 1]
				
				if(np.sum(dat_hB0) > 1): #if actual data exists
					plt.imshow(dat_hB0, cmap = 'gray', aspect = 'auto',\
							vmin = vmin, vmax = vmax, \
							extent = (wg[idx_observed - rang], \
								wg[idx_observed + rang], dat.shape[0], 0))
					plt.grid(False)
			
			elif(((hB0 + crdelt1*rang)*(1 + z[redshift_idx]) > wg[-1]) & \
			((hB0 + cdelt_red*rang)*(1 + z[redshift_idx]) < wg_red[-1])):
				idx_observed = np.abs(wg_red - observed_wavelength).argmin()
				dat_hB0 = dat_red[:, idx_observed - rang: idx_observed + rang + 1]
				
				if(np.sum(dat_hB0) > 1): #if actual data exists
					plt.imshow(dat_hB0, cmap = 'gray', aspect = 'auto',\
						vmin = vmin, vmax = vmax, \
						extent = (wg_red[idx_observed - rang], \
								wg_red[idx_observed + rang], dat_red.shape[0], 0))
					plt.grid(False)
		else: #if redder mask
			if((hB0 + crdelt1*rang)*(1 + z[redshift_idx]) < wg[-1]): #make sure rightmost point exists within the wg
				idx_observed = np.abs(wg - observed_wavelength).argmin()
				dat_hB0 = dat[:, idx_observed - rang: idx_observed + rang + 1]
				
				if(np.sum(dat_hB0) > 1): #if actual data exists
					plt.imshow(dat_hB0, cmap = 'gray', aspect = 'auto',\
							vmin = vmin, vmax = vmax, \
							extent = (wg[idx_observed - rang], \
								wg[idx_observed + rang], dat.shape[0], 0))
					plt.grid(False)
		
		# OIII1
		observed_wavelength = OIII1*(1 + z[redshift_idx])
		
		plt.subplot2grid((4,4), (0,1))
		plt.title('[O III] 49')
		if(maskred_flag): #if bluer mask
			if((OIII1 + crdelt1*rang)*(1 + z[redshift_idx]) < wg[-1]): #make sure rightmost point exists within the wg
				idx_observed = np.abs(wg - observed_wavelength).argmin()
				dat_OIII1 = dat[:, idx_observed - rang: idx_observed + rang + 1]
				
				if(np.sum(dat_OIII1) > 1): #if actual data exists
					plt.imshow(dat_OIII1, cmap = 'gray', aspect = 'auto',\
							vmin = vmin, vmax = vmax, \
							extent = (wg[idx_observed - rang], \
								wg[idx_observed + rang], dat.shape[0], 0))
					plt.grid(False)
			
			elif(((OIII1 + crdelt1*rang)*(1 + z[redshift_idx]) > wg[-1]) & \
			((OIII1 + cdelt_red*rang)*(1 + z[redshift_idx]) < wg_red[-1])):
				idx_observed = np.abs(wg_red - observed_wavelength).argmin()
				dat_OIII1 = dat_red[:, idx_observed - rang: idx_observed + rang + 1]
				
				if(np.sum(dat_OIII1) > 1): #if actual data exists
					plt.imshow(dat_OIII1, cmap = 'gray', aspect = 'auto',\
						vmin = vmin, vmax = vmax, \
						extent = (wg_red[idx_observed - rang], \
								wg_red[idx_observed + rang], dat_red.shape[0], 0))
					plt.grid(False)
		else: #if redder mask
			if((OIII1 + crdelt1*rang)*(1 + z[redshift_idx]) < wg[-1]): #make sure rightmost point exists within the wg
				idx_observed = np.abs(wg - observed_wavelength).argmin()
				dat_OIII1 = dat[:, idx_observed - rang: idx_observed + rang + 1]
				
				if(np.sum(dat_OIII1) > 1): #if actual data exists
					plt.imshow(dat_OIII1, cmap = 'gray', aspect = 'auto',\
							vmin = vmin, vmax = vmax, \
							extent = (wg[idx_observed - rang], \
								wg[idx_observed + rang], dat.shape[0], 0))
					plt.grid(False)
		
		#OIII2
		observed_wavelength = OIII2*(1 + z[redshift_idx])
		
		plt.subplot2grid((4,4), (0,2))
		plt.title('[O III] 50')
		if(maskred_flag): #if bluer mask	
			if((OIII2 + crdelt1*rang)*(1 + z[redshift_idx]) < wg[-1]): #make sure rightmost point exists within the wg
				idx_observed = np.abs(wg - observed_wavelength).argmin()
				dat_OIII2 = dat[:, idx_observed - rang: idx_observed + rang + 1]
				
				if(np.sum(dat_OIII2) > 1): #if actual data exists
					plt.imshow(dat_OIII2, cmap = 'gray', aspect = 'auto',\
							vmin = vmin, vmax = vmax, \
							extent = (wg[idx_observed - rang], \
								wg[idx_observed + rang], dat.shape[0], 0))
					plt.grid(False)
			
			elif(((OIII2 + crdelt1*rang)*(1 + z[redshift_idx]) > wg[-1]) & \
			((OIII2 + cdelt_red*rang)*(1 + z[redshift_idx]) < wg_red[-1])):
				idx_observed = np.abs(wg_red - observed_wavelength).argmin()
				dat_OIII2 = dat_red[:, idx_observed - rang: idx_observed + rang + 1]
				
				if(np.sum(dat_OIII2) > 1): #if actual data exists
					plt.imshow(dat_OIII2, cmap = 'gray', aspect = 'auto',\
							vmin = vmin, vmax = vmax, \
							extent = (wg_red[idx_observed - rang], \
									wg_red[idx_observed + rang], dat_red.shape[0], 0))
					plt.grid(False)
		else: #if redder mask
			if((OIII2 + crdelt1*rang)*(1 + z[redshift_idx]) < wg[-1]): #make sure rightmost point exists within the wg
				idx_observed = np.abs(wg - observed_wavelength).argmin()
				dat_OIII2 = dat[:, idx_observed - rang: idx_observed + rang + 1]
				
				if(np.sum(dat_OIII2) > 1): #if actual data exists
					plt.imshow(dat_OIII2, cmap = 'gray', aspect = 'auto',\
							vmin = vmin, vmax = vmax, \
							extent = (wg[idx_observed - rang], \
								wg[idx_observed + rang], dat.shape[0], 0))
					plt.grid(False)
		
		#spectra
		plt.subplot2grid((4,4), (1,0), colspan=4, rowspan=1)
		plt.ylabel('spectra')
		plt.plot(wg_windowed, image_windowed)
		plt.plot(wg_windowed, model, c = 'k', linestyle = "--")
		plt.ylim([np.min(model) - 1., np.max(model) + 1]) # ignore outliers from ruining plot scale
		
		#del chi2
		plt.subplot2grid((4,4), (2,0), colspan=4, rowspan=1)
		plt.ylabel(r'$\Delta \chi^2$')
		plt.plot(lambda0*(1 + z_range), delChi2data[slit_idx][w_idx])
		
		#ivar
		plt.subplot2grid((4,4), (3,0), colspan=4, rowspan=1)
		plt.ylabel('ivar')
		plt.plot(wg_windowed, ivar_windowed)
		plt.ylim(ymin = 0)
		
		# title
		SNRmax = SNRdata[idx, w_idx, redshift_idx]
		SNRmax = np.round(SNRmax, 2)	
		plt.suptitle(f'Mask: {maskname} Slit: {slit_idx} z : {np.round(z[redshift_idx],2)} SNR: {SNRmax} vel: {sigma_v[w_idx]} km/s rel_str: {np.round(rel_str,2)}', \
						fontsize = 15, fontname = 'serif')
		
		# fit subplots and save fig
		fig.tight_layout()
		fig.set_size_inches(w=11,h=7)
		fig_name = '../results/PeakZoom/' + str(maskname) + '/' + str(maskname) + "-Slit-" \
		+ str(slit_idx) + '.png'
		fig.savefig(fig_name, dpi = 250, bbox_inches = 'tight')
		plt.close()

#user input of maskname; as str
maskname = sys.argv[1]
flag = (sys.argv[2].lower() == 'true') #user entered input whether bluer or redder mask

#read in data based on maskname
data = datareader(maskname)
image = data['data_ivar'][:, 0, :]
ivar = data['data_ivar'][:, 1, :]
datarows = len(image)

zdata, widthsdata, SNRdata, Ampdata, delChi2data = npy_reader(maskname)

#best models indices from best-model.py
bestmodels = np.load("../results/Max_z_n_width/" + maskname + "-indices.npy")

_, wg = wave_grid(data)

from time import time

start = time()
for idx in range(1, datarows):
	redshift_idx = bestmodels[idx-1, 0]
	w_idx = bestmodels[idx-1, 1]
	
	if(redshift_idx != -999): #ignore nans
		if(SNRdata[idx, w_idx, redshift_idx] >= 10):
			PeakZoom(maskname = maskname, slit_idx = idx, z = zdata, widths = widthsdata, \
			SNRdata = SNRdata, Ampdata = Ampdata, delChi2data = delChi2data, \
			image = image, ivar = ivar, wg = wg, redshift_idx = redshift_idx, \
			w_idx = w_idx, rel_str = relative_strength, fudge_fac = fudge_factor, \
			maskred_flag = flag)
		 
end = time()
print(f'Time: {end - start}')