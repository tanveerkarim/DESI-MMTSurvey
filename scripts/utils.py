"""This module contains all the necessary functions to do MMT BinoSpec 1D analysis.
Authors: Tanveer Karim and Jae Lee
Latest version: 09-Jul-2018"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.ndimage.filters import median_filter
from global_var import *
import pandas as pd
import os
plt.style.use('ggplot')

"""------START OF TANVEER'S CODE------"""

def datareader(maskname, dir_name = "../npz_files/"):
	"""Reads mask data for use by the other functions in this module
	Parameters
	----------
	maskname: name of the mask + '-' + grating number
	"""
	
	fname = maskname + '-' + 'spec1d.npz'
	data = np.load(dir_name + fname)
	
	return data

def npy_reader(maskname):
	"""reads in all the .npy products from outputdata folder"""
	
	if(not isinstance(maskname, str)):
		maskname = str(maskname)
		
	zdata = np.load("../results/outputdata/" + maskname + "/" + maskname + "-0.7-z_range.npy")
	widthsdata = np.load("../results/outputdata/" + maskname + "/" + maskname + "-0.7-widths.npy")
	SNRdata = np.load("../results/outputdata/" + maskname + "/" + maskname + "-0.7-SNR.npy")
	Ampdata = np.load("../results/outputdata/" + maskname + "/" + maskname + "-0.7-Amp.npy")
	delChi2data = np.load("../results/outputdata/" + maskname + "/" + maskname + "-0.7-delChi2.npy")
	
	return zdata, widthsdata, SNRdata, Ampdata, delChi2data

def wave_grid(data):
	"""Returns wavegrid based on header file from data"""
	
	crval1 = float(str(data['headers'][1]).split("CRVAL1")[1].\
				split("=")[1].split("/")[0]) #Starting value
	cdelt1 = float(str(data['headers'][1]).split("CDELT1")[1].\
				split("=")[1].split("/")[0]) #Pixel size
	
	collapsedSpectrum = data['data_ivar'][:, 0, :]
	
	wave_grid = crval1 + cdelt1 * np.arange(collapsedSpectrum[1].shape[0])
	wave_grid *= 10 #Convert wave_grid to Angstrom from nm
	crval1 *= 10
	
	return crval1, wave_grid

def lambda_to_z(wavelength):
	"""Converts wavelength grid to redshift grid"""
	
	return (wavelength/lambda0 - 1)

def Window(z, wg, window_width = 0.005):
	"""Returns a range of pixel in the specified window width
	
	Parameters
	----------
	z: Centre of the window
	wg: wave grid that needs to be windowed
	z_grid: redshift grid of the wave_grid
	window_width: size of the window in redshift space
	
	Returns
	-------
	windowed_array: windowed array of the windowing_array	
	"""

	zgrid = lambda_to_z(wg)
	
	#Multiply window_width by (1+z) to maintain same window size in lambda space
	windowed_array = wg[(zgrid > (z - window_width*(1+z))) \
	& (zgrid < (z + window_width*(1+z)))]
	
	return windowed_array

def Model(z, wg, width, relative_strength, Amp = 1.):
	"""Returns Gaussian filter models at redshift z for
	a number of different Gaussian widths
	
	Parameters
	----------
	z: redshift at which the model is being tested
	wg: pixel grid of the Window
	width: width array of the Gaussian doublets
	relative_strength: Relative strength of 3729 wrt 3727
	Amp: amplitude of the Gaussian doublets
	
	Returns
	--------
	model: Gaussian models in the range of [z - window_width, z + window_width]
	"""
	
	"""broadcasting x because both wg2 and width are 1d arrays so in order to get a model 
	of the form (len(wg2), len(width)), need to broadcast wg2"""
	Gaussian = lambda x, mean, std: (1/np.sqrt(2*np.pi*std**2))*\
	np.exp(-0.5*((x[:, np.newaxis] - mean)/std)**2)
	#relative_strength = 1. #http://www.ucolick.org/~simard/phd/root/node21.html
	
	model = (Amp/(1+relative_strength))*(relative_strength*Gaussian(wg, lambda27*(1+z), width)\
	+ Gaussian(wg, lambda29*(1+z), width))
		
	return model

def SNR_calculator(maskname, data, rel_strngth, fudge_factor):

	"""maskname[-3:] yields the grating number. z_range changes depending
	on maskname because of variation in grating. The start and end points
	are chosen by inspecting the header file of the data."""
	
	#Read data
	image = data['data_ivar'][:, 0, :]
	ivar = data['data_ivar'][:, 1, :]
	crval1, wg = wave_grid(data)
	z_grid = lambda_to_z(wg) #Convert wavelength space to redshift space
	
	#set max value for z_range
	zmaximum = (10000./lambda0) - 1. #max lambda is 10000 A
	zmaximum_idx = np.abs(zmaximum - z_grid).argmin()
	zmaximum = z_grid[zmaximum_idx]
	
	#Set redshift hypothesis grid
	z_range = np.arange(z_grid[0], zmaximum, 0.0002)
	
	"""Gaussian width, sigma = sqrt(sigma_lambda^2 + sigma_slit^2) where, 
	sigma_lambda = sigma_v/c*lambda(z); sigma_v = [0, 300] km/s
	sigma_slit = 3.3/sqrt(12)*delLambda_pixel	
	"""
	
	delLambda_pixel = float(str(data['headers'][1]).split("CDELT1")[1]\
		.split("=")[1].split("/")[0])*10. #size of the pixel in angstrom
	
	#fudge factor introduced to make fit better. Required because size of object may be smaller
	sigma_slit = ((3.3/np.sqrt(12))*delLambda_pixel)*fudge_factor 
	c = 299792.458 #km/s
			
	def widthlist(z):
		"""Returns an array of possible Gaussian widths for the [OII] doublet
		model testing"""
		
		sigma_lambda = sigma_v/c*(lambda0*(1 + z)) #in observing frame
	
		return np.sqrt(sigma_lambda**2 + sigma_slit**2)
	
	#sigma_v size same as number of Gaussian width models
	#store all the SNR values (z x num. of objects x width)
	SNRs = np.zeros((z_range.size, image.shape[0], sigma_v.size))
	
	#Save all the amplitudes and chi_sq to pass this to the PeakZoom function
	#Amps -> (z x num. of objects x width)
	#del_chi_sq -> (z x num. of objects x width)
	Amps = np.zeros(SNRs.shape) 
	del_chi_sq = np.zeros(SNRs.shape)
	
	#Store all the widths b/c widths are a func. of z
	widths = np.zeros((z_range.size, sigma_v.size))
	
	#-------------------------------------------#
	#WORK ON THIS LATER
	#Save all the medians to pass to PeakZoom function
	medians = np.zeros((z_range.size, image.shape[0]))
	#-------------------------------------------#
	
	#for i, z in tqdm(enumerate(z_range)):
	for i, z in enumerate(z_range):
		#print(i)
		wg2 = Window(z, wg)
		widths[i] = widthlist(z) #Annoying bug. MUST save widths becuase widths = widths(z). Previously was using last z 
		model = Model(z, wg2, widths[i], relative_strength = rel_strngth)		
		
		#Find the idx of the edges of the windows and slice the image file to multiply with modelPrime
		minidx = np.where(wg == np.min(wg2))[0][0] 
		maxidx = np.where(wg == np.max(wg2))[0][0]
		imageSliced = np.copy(image[:,minidx:maxidx+1]) #imageSliced -> (num. of obj x window wavelength)
		ivarSliced = np.copy(ivar[:,minidx:maxidx+1]) #ivarSliced -> (num. of obj x window wavelength)
		
		#create a copy of imageSliced to set peak area to 0 to get proper dc level
		#[z-0.001, z + 0.001] range where we expect the doublet to be and set it to 0
		upperlam = lambda0*(1 + (z + 0.001)) 
		lowerlam = lambda0*(1 + (z - 0.001))
		
		#find corresponding indices to set imageSliced values inside that range to be 0
		rightidx = np.abs(wg2 - upperlam).argmin()
		leftidx = np.abs(wg2 - lowerlam).argmin()
		
		imageSlicedtmp = np.copy(imageSliced)
		imageSlicedtmp[:, leftidx:rightidx+1] = 0
		
		#Ignore 0s in imageSliced when calculating median 
		#Source: https://stackoverflow.com/questions/22049140/how-can-i-ignore-zeros-when-i-take-the-median-on-columns-of-an-array/22049849#22049849
		median_val = np.apply_along_axis(lambda v: np.median(v[v!=0]), 1, imageSlicedtmp)
		median_val[np.isnan(median_val)]=0. #Need to do this because some windows are all 0 and results in NaN in np.median(v[v!=0])
		medians[i] = median_val
		imageSlicedtmp2 = imageSliced - median_val[:, np.newaxis] #median subtraction
		imageSliced = imageSlicedtmp2 
		
		"""numpy dot -> sum product over last axis of a and
		second-to-last axis of b where np.dot(a,b)
		Here, (imageSliced*ivarSliced) -> elementwise multiplication
		producing (num of obj x window) and model -> (window x width).
		Hence, summing over window, i.e. range of pixels gives us the 
		numerator -> (num of obj x width)"""
		trm1 = np.multiply(imageSliced, ivarSliced)[:,:, np.newaxis]
		trm2 = model[np.newaxis, :, :]
		Numerator = trm1 * trm2
		Numerator = np.sum(Numerator, axis = 1)
		
		"""numpy dot -> sum product over last axis of a and
		second-to-last axis of b where np.dot(a,b)
		Here, (model * model) -> elementwise multiplication
		producing (window x width) and ivar -> (num of obj x window).
		Hence, summing over window, i.e. range of pixels gives us the 
		denominator. We can do ivar.model^2 because since the operation
		is pixel-wise, it does not matter whether it is model^2.ivar
		or ivar.model^2.
		Denominator -> (num of obj x width)"""
		trm3 = ivarSliced[:, :, np.newaxis]
		trm4 = np.multiply(model, model)
		trm5 = trm4[np.newaxis, :, :]
		Denominator = trm3 * trm5
		Denominator = np.sum(Denominator, axis = 1)
		
		"""
		M' = M/sigma_px; D' = D/sigma_px
		A = (D'.M')/(M'.M')
		sigmaA^(-2) = M'.M'
		Let, isigmaA = sqrt(M'.M')
		SNR = A/sigmaA => SNR = A*isigmaA
		"""
		
		Amp = Numerator/(Denominator+1.0e-100) #
		isigmaA = np.sqrt(Denominator)
		SNR = Amp*isigmaA
	
		SNRs[i] = SNR
		Amps[i] = Amp
		
		"""chi_sq = sum over pixel ((image - model*amp)^2*ivar)
		del chi_sq = sum over pixel ((image - model*amp)^2*ivar) - (image^2*ivar)
		firstrm = image
		secondtrm = model*amp
		thirdtrm = image^2 = firstrm^2
		delterm1 = ((image - model*amp)^2*ivar)
		delterm2 = (image^2*ivar)
		"""
		
		firstrm = imageSliced[:, np.newaxis, :] #(ngal x 1 x range)
		secondtrm = Amp[:,:,np.newaxis] * model.T[np.newaxis, :, :] #(ngal x width x range)
		thirdtrm = firstrm**2 #(ngal x 1 x range)
		ivartrm = ivarSliced[:, np.newaxis, :] #(ngal x 1 x range)
		
		diff = firstrm - secondtrm
		delterm1 = (diff**2)*ivartrm
		delterm2 = thirdtrm*ivartrm
		
		del_chi = delterm1 - delterm2
		del_chi2 = np.nansum(del_chi, axis = 2)
		del_chi_sq[i] = del_chi2
		
	SNRs_final = SNRs.transpose([1, 2, 0]) #This maintains the indices
	Amps_final = Amps.transpose([1, 2, 0])
	del_chi_sq_final = del_chi_sq.transpose([1, 2, 0])
	
	return z_range, widths, SNRs_final, Amps_final, del_chi_sq_final, medians 

def plotterSpectra1D(maskname, data, idx):
	"""Returns plots of the 1d spectra and the inverse variance
	
	Parameters
	----------
	idx: Index number of the slit; ignore idx = 0
	"""
	
	image = data['data_ivar'][:, 0, :]
	ivar = data['data_ivar'][:, 1, :]
	
	imagetmp = image[idx, :]
	ivartmp = ivar[idx, :]
	
	#find the edges of the 1d spectra
	"""image_copy = pd.DataFrame(imagetmp)
	idx_min = image_copy.first_valid_index()
	idx_max = image_copy.last_valid_index()
	if(idx_min > 5):
		idx_min = idx_min + 5
	if(idx_max < len(image_copy) - 6):
		idx_max = idx_max + 5"""
	
	#wavegrid
	wavegrid = wave_grid(data)
	
	f, axarr = plt.subplots(2, sharex=True)
	axarr[0].plot(wavegrid, imagetmp)
	#axarr[0].set_xlim(wavegrid[idx_min], wavegrid[idx_max]) #set limit to show only good data
	axarr[0].set_ylim(-2, 3)
	axarr[0].set_title('Mask: ' + maskname + ', ' + 'Slit ' + str(idx) + "\n" + "1D spectra" \
					,  fontsize = 15, fontname = 'serif')
	axarr[1].plot(wavegrid, ivartmp)
	#axarr[1].set_xlim(wavegrid[idx_min], wavegrid[idx_max]) #set limit to show only good data
	axarr[1].set_title('1D inverse variance', fontsize = 10, fontname = 'serif')
	plt.savefig('results/spectra1d/' + maskname + '/' + maskname + '-' + str(idx) + '-spectra1d.pdf', dpi = 250, bbox_inches = None)
	plt.close()

def plotter2D(maskname, idx, z, widths, SNRdata, delChi2data):
	"""Returns redshift vs. width plot with SNR strength and del_chi2 strength.
	
	Parameters
	----------
	maskname: name of the mask + '-' + grating number
	idx: index of a slit for a given maskname
	z: 0th output of the SNR_calculator function; redshift range
	widths: 1st output of the SNR_calculator function; width range
	SNRdata: 2nd output of the SNR_calculator function; SNR data cube
	delChi2data: 3rd output of the SNR calculator function; delChi2 data cube
	
	Returns
	-------
	PDF of the SNR and delChi2 2D plots
	"""
	
	plt.imshow(SNRdata[idx], aspect='auto', interpolation='None', \
			extent=[np.min(z), np.max(z), np.min(widths), np.max(widths)], vmin=0)#, vmax=7)
	plt.colorbar()
	plt.ylabel('width', fontsize = 15, fontname = 'serif')
	plt.xlabel('redshift', fontsize = 15, fontname = 'serif')
	plt.title('Mask: ' + maskname + ', ' + 'Slit ' + str(idx),  fontsize = 15, fontname = 'serif')
	plt.savefig("../results/SNR2D/"  + maskname + '/' + maskname + '-' + str(idx) + \
	"-SNR2d.pdf", dpi = 250, bbox_inches = None)
	plt.close()
	
	#Plot del_chi2
	plt.imshow(delChi2data[idx], aspect='auto', interpolation='None', \
			extent=[np.min(z), np.max(z), np.min(widths), np.max(widths)])
	plt.colorbar()
	plt.ylabel('width', fontsize = 15, fontname = 'serif')
	plt.xlabel('redshift', fontsize = 15, fontname = 'serif')
	plt.title('Mask: ' + maskname + ', ' + 'Slit ' + str(idx),  fontsize = 15, fontname = 'serif')
	plt.savefig("../results/delChi22D/"  + maskname + '/' + maskname + '-' + str(idx) + \
	"-delChi22d.pdf", dpi = 250, bbox_inches = None)
	plt.close()

def SNRvz(maskname, idx, z, widths, SNRdata, Ampdata, delChi2data, medians, image, ivar, wg, rel_str, fudge_fac, maskred_flag):
	"""Returns SNR vs z plot per slit and redshift and w values
	Parameters
	----------
	maskname: name of the mask + '-' + grating number
	idx: index of a slit for a given maskname
	z: 0th output of the SNR_calculator function; redshift range
	widths: 1st output of the SNR_calculator function; width range
	SNRdata: 2nd output of the SNR_calculator function; SNR data cube
	delChi2data: 3rd output of the SNR calculator function; delChi2 data cube
	image: spectra 1d -> pass to PeakZoom
	ivar: inverse variance 1d -> pass to PeakZoom
	wg: wavelength grid of maskname
	rel_str: relative strength of [O II] 27 wrt 29
	fudge_fac: fudge factor for sigma_slit
	maskred_flag: whether a red mask of maskname exists; if so, check that 
					for plotting purposes
	"""	
	

	#Find width and z indices for highest SNR
	w_idx, redshift_idx = np.argwhere(delChi2data[idx] == np.min(delChi2data[idx]))[0]
	
	if(SNRdata[idx, w_idx, redshift_idx] >= 10):
		#Generate SNR vs z plot per slit
		#plt.plot(z, SNRdata[idx, w_idx])
		#plt.axhline(7, c = 'red') #Threshold of SNR = 7
		#plt.ylabel('SNR', fontsize = 15, fontname = 'serif')
		#plt.xlabel('redshift', fontsize = 15, fontname = 'serif')
		#plt.title('Mask: ' + maskname + ', ' + 'Slit ' + str(idx) +"\n" +\
		#	"z = " + str(np.round(z[redshift_idx], 3)) + ', w = '\
		#	+ str(np.round(widths[w_idx],2)), fontsize = 15, fontname = 'serif')
		#plt.xlim([z[redshift] - .1, z[redshift] + .1])
		#plt.savefig("../results/SNRvsRedshift/"  + maskname + '/' + maskname\
		#+ '-' + str(idx) + "-SNR_vs_z.png", dpi = 250, bbox_inches = None)
		#plt.close()
		
		#PeakZoom(int(maskname), idx, z, widths, SNRdata, Ampdata, delChi2data,\
		#medians, image, ivar, wg, redshift_idx, w_idx,\
		#rel_str = rel_str, fudge_fac = fudge_fac, maskred_flag = maskred_flag)
		hypothesis_tester(maskname, idx, z, widths, \
			image, ivar, wg, redshift_idx, w_idx, rel_str, fudge_fac, maskred_flag)
		return z[redshift_idx], sigma_v[w_idx]
	else:
		return np.nan, np.nan

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
	
	#Find width and z indices for highest neg delta chi2 change
	w_idx, redshift_idx = np.argwhere(delChi2data[idx] == np.min(delChi2data[idx]))[0]
	
	if(SNRdata[idx, w_idx, redshift_idx] >= 10):
		return z[redshift_idx], sigma_v[w_idx], redshift_idx, w_idx
	else:
		return -999, -999, -999, -999
		
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
		hB0 = 4862.68
		OIII1 = 4960.295
		OIII2 = 5008.240   
		
		crdelt1 = 0.62 #ad hoc value crdelt1*rang because rang is in px 
            #space so convert px to wavelength space and crdel1 is the 
            #conversion factor
    
		"""set range that will plot same window width as that of [O II] doublet.
		given the index of a given line, we will generate 2D plots of 
		xlim = [idx - crdelt1*rang, idx + crdelt1*rang]"""
		rang = np.median(np.arange(0, wg_windowed.shape[0]))
		rang = np.int(rang)
		
		#Plot figure with subplots of different sizes
		alpha = 0.7 #transparency of vertical lines
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

def hypothesis_tester(maskname, slit_idx, z, widths, \
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
		data_err, list_headers = preprocess_bino(masknumber=maskname)
		dat, err, header = extract_single_data(data_err, list_headers, idx)
		if(maskred_flag == True):
			redmask = np.str(np.int(maskname) + 1) #red mask is one up from blue mask
			data_err, list_headers = preprocess_bino(masknumber=redmask)
			dat_red, _, _ = extract_single_data(data_err, list_headers, idx)
			
			tmpdata = datareader(redmask)
			crval_red, wg_red = wave_grid(tmpdata) 
			cdelt_red = tmpdata['headers'][1]['CDELT1']*10 #convert to angstrom
		else: #if redder mask is default
			bluemask = np.str(np.int(maskname) - 1) #blue mask is one down from red mask
			data_err, list_headers = preprocess_bino(masknumber=bluemask)
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
		#Values from SDSS
		OII = lambda0
		hB0 = 4862.68
		OIII1 = 4960.295
		OIII2 = 5008.240   
		ha = 6564.61
		
		rang = np.median(np.arange(0, wg_windowed.shape[0]))
		rang = np.int(rang) # rang is in px space
		
		#Plot figure with subplots of different sizes
		alpha = 0.7
		crdelt1 = 0.62 #ad hoc value crdelt1*rang because rang is in px 
				#space so convert px to wavelength space and crdel1 is the 
				#conversion factor
		vmin = -0.5
		vmax = 1
		
		fig = plt.figure(1, figsize=(15,15))
		# set up subplot grid
		subgrid_plot_size = (5,5)
		gridspec.GridSpec(5,5)
		aspect = 'auto'
		fontsize = 25
		
		def plot_hypothesis(emission_line, z_hypothesis, grid_loc, maskred_flag = maskred_flag):
			species = dict({hB0:r'H$\beta$', OIII1:'[OIII] 49',
						OIII2:'[OIII] 50', ha:r'H$\alpha$'})
				
			observed_wavelength = emission_line*(1 + z_hypothesis)
			
			if(maskred_flag): #if testing for bluer mask
				if((emission_line + crdelt1*rang)*(1 + z_hypothesis) < wg[-1]): #make sure rightmost point exists within the wg
					idx_observed = np.abs(wg - observed_wavelength).argmin()
					dat_emission_line = dat[:, idx_observed - rang: idx_observed + rang + 1]
				
					if(np.sum(dat_emission_line) > 1): #if actual data exists
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
					
				elif(((emission_line + crdelt1*rang)*(1 + z_hypothesis) > wg[-1]) & \
					((emission_line + cdelt_red*rang)*(1 + z_hypothesis) < wg_red[-1])):
					idx_observed = np.abs(wg_red - observed_wavelength).argmin()
					dat_emission_line = dat_red[:, idx_observed - rang: idx_observed + rang + 1]
				
					if(np.sum(dat_emission_line) > 1): #if actual data exists
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
						
				else:
					if(grid_loc[0] == 0):
							plt.subplot2grid(subgrid_plot_size, grid_loc)
							plt.axis('off')
							plt.title(species[emission_line], fontsize = fontsize)
							
			else: #for red masks    
				if(((emission_line + crdelt1*rang)*(1 + z_hypothesis) < wg[0]) & \
					((emission_line + cdelt_blue*rang)*(1 + z_hypothesis) > wg_blue[0])): #check bluer mask
					idx_observed = np.abs(wg_blue - observed_wavelength).argmin()
					dat_emission_line = dat_blue[:, idx_observed - rang: idx_observed + rang + 1]
				
					if(np.sum(dat_emission_line) > 1): #if actual data exists
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
				
				elif((emission_line + crdelt1*rang)*(1 + z_hypothesis) < wg[-1]): #make sure rightmost point exists within the wg
					idx_observed = np.abs(wg - observed_wavelength).argmin()
					dat_emission_line = dat[:, idx_observed - rang: idx_observed + rang + 1]
				
					if(np.sum(dat_emission_line) > 1): #if actual data exists
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
		#plt.xticks([])
		
		# OIII1
		#plt.subplot2grid(subgrid_plot_size, (0,2))
		plt.axis('off')
		plot_hypothesis(OIII1, z[redshift_idx], (0,2))
		#plt.title('[O III] 49')
		#plt.xticks([])
		plt.yticks([])
		
		#OIII2
		#plt.subplot2grid(subgrid_plot_size, (0,3))
		plot_hypothesis(OIII2, z[redshift_idx], (0,3))
		#plt.title('[O III] 50')
		#plt.xticks([])
		plt.yticks([])
		
		# hB
		#plt.subplot2grid(subgrid_plot_size, (0,1))
		plot_hypothesis(hB0, z[redshift_idx], (0,1))
		#plt.title(r'H$\beta$')
		#plt.xticks([])
		plt.yticks([])
		
		# ha
		#plt.subplot2grid(subgrid_plot_size, (0,4))
		plot_hypothesis(ha, z[redshift_idx], (0,4))
		#plt.title(r'H$\alpha$')
		#plt.xticks([])
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
		#plt.xticks([])
		plt.yticks([])
		plt.grid(False)
		
		# OIII1
		#plt.subplot2grid(subgrid_plot_size, (1,2))        
		plot_hypothesis(OIII1, z_hypothesis, (1,2))
		#plt.xticks([])
		plt.yticks([])
			
		#OIII2
		#plt.subplot2grid(subgrid_plot_size, (1,3))
		plot_hypothesis(OIII2, z_hypothesis, (1,3))
		#plt.xticks([])
		plt.yticks([])
		
		# OII
		#plt.subplot2grid(subgrid_plot_size, (1,0))
		plot_hypothesis(OII, z_hypothesis, (1,0))
		#plt.xticks([])
		
		# ha
		#plt.subplot2grid(subgrid_plot_size, (1,4))
		plot_hypothesis(ha, z_hypothesis, (1,4))
		#plt.xticks([])
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
		#plt.xticks([])
		plt.yticks([])
		plt.grid(False)
		
		# Hbeta
		#plt.subplot2grid(subgrid_plot_size, (2,1))        
		plot_hypothesis(hB0, z_hypothesis, (2,1))
		#plt.xticks([])
		plt.yticks([])
			
		#OIII2
		#plt.subplot2grid(subgrid_plot_size, (2,3))
		#plt.xticks([])
		plt.yticks([])
		plot_hypothesis(OIII2, z_hypothesis, (2,3))
		
		# OII
		#plt.subplot2grid(subgrid_plot_size, (2,0))
		plot_hypothesis(OII, z_hypothesis, (2,0))
		#plt.xticks([])
		
		# ha
		#plt.subplot2grid(subgrid_plot_size, (2,4))
		plot_hypothesis(ha, z_hypothesis, (2,4))
		#plt.xticks([])
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
		#plt.subplot2grid(subgrid_plot_size, (3,1))        
		plot_hypothesis(hB0, z_hypothesis, (3,1))
		plt.yticks([])
			
		#OIII1
		#plt.subplot2grid(subgrid_plot_size, (3,2))
		plot_hypothesis(OIII1, z_hypothesis, (3,2))
		plt.yticks([])
		
		# OII
		#plt.subplot2grid(subgrid_plot_size, (3,0))
		plot_hypothesis(OII, z_hypothesis, (3,0))
		
		# ha
		#plt.subplot2grid(subgrid_plot_size, (3,4))
		plot_hypothesis(ha, z_hypothesis, (3,4))
		#plt.xticks([])
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
		#plt.subplot2grid(subgrid_plot_size, (4,1))        
		plot_hypothesis(hB0, z_hypothesis, (4,1))
		plt.yticks([])
			
		#OIII1
		#plt.subplot2grid(subgrid_plot_size, (4,2))
		plot_hypothesis(OIII1, z_hypothesis, (4,2))
		plt.yticks([])
		
		# OII
		#plt.subplot2grid(subgrid_plot_size, (4,0))
		plot_hypothesis(OII, z_hypothesis, (4,0))
		
		# OIII2
		#plt.subplot2grid(subgrid_plot_size, (4,3))
		plot_hypothesis(OIII2, z_hypothesis, (4,3))
		#plt.xticks([])
		plt.yticks([])
		
		
		plt.suptitle(f'Mask: {maskname} Slit: {slit_idx} z : {np.round(z[redshift_idx],2)} vel: {sigma_v[w_idx]} km/s', \
						fontsize = fontsize, fontname = 'serif')
		plt.tight_layout(rect=[0, 0.03, 1, 0.95])
		plt.subplots_adjust(wspace=0.1, hspace=0.1)
		
		#plt.tight_layout(rect=[0, 0.03, 1, 0.95]) #so that main title and subtitles don't overlap
		#plt.subplots_adjust(wspace=0.05, hspace=0.05)
		#fig.set_size_inches(w=11,h=7)
		fig_name = '../results/hyp_test/' + str(maskname) + '/' + str(maskname) + "-Slit-" \
		+ str(slit_idx) + '.png'
		fig.savefig(fig_name, dpi = 250, bbox_inches = 'tight')
		plt.close()
	
def idx_peakOII(wavegrid, redz, idx_min=0, idx_max=None):
    """
    Given a wavelength grid and a redshift, return the indices corresponding to
    the 3727 [O II] peak
    """

    OII = 3727
    peak_list = [OII]
    
    if idx_max is None:
        idx_max = wavegrid.size-1
    
    # Compute redshifted location
    peak_redshifted = OII*(1+redz)    
        
    # Compute wavegrid index corresponding to the location. Return -1 if outside the bound.
    idx = find_nearest_idx(wavegrid, peak_redshifted)
    if (idx >=idx_min) and (idx < idx_max):
            index = idx
    else:
        index = -1

    return index

"""------START OF JAE'S CODE------"""

def preprocess_bino(data_dir, \
	fname_data = "obj_abs_slits_lin.fits", fname_err = "obj_abs_err_slits_lin.fits",\
					fname_extr = "obj_abs_slits_extr.fits", masknumber = 1624, pos_width = 10): 
	"""
	Based on Jae's original code. Goal is to process binospec data in a way
	such that it can be used to extract 1D spectra. 
	
	1. For a given mask, finds where the galaxies are along the position
	axis. This is required because we extended the slits whenever possible.
	OBPOSVAL from the header tells us where the central points of the 
	galaxies per slit are so that we can do extraction only around that 
	region.
	
	2. Extract +/- pos_width pixels around the galaxies for 1D extraction. 
	Default set to 10 px.
	
	3. The native data unit is ergs/cm^2/s/nm. Preprocessor changes this to
	10^-17 ergs/cm^2/s/Angstrom.
	
	4. Following conditions imposed -
	Data: 
		- If NaN: Data ZERO and Error infinity.
	Error:
		- If NaN: Error infinity and Data ZERO.
		
	Args:
		fname_data: data file name
		fname_err: error file name
		fname_extr: head file name to get OBPOSVAL
		masknumber: mask number
		post_width: amount of pixels to be used for extraction above and 
			below a given galaxy
	Returns: 
		data_err: A numpy array of shape (Ntargets+1, 2, pos_width, num_cols). 
				- Though the native data have different number of rows, we use a single fixed number here.
				- Channel 0: Data
				- Channel 1: Error
		list_headers: List of headers
	
	Note, first spectrum in native data is saved in loc "1". We follow the same convention.
	"""
	
	data = fits.open(data_dir + fname_data)
	err = fits.open(data_dir + fname_err)
	extr = fits.open(data_dir + fname_extr)
	extr = extr[1] #Need only one header to get all the OBPOSVALs
	
	infinity = 1e60
	#unit_conversion = 10**19
	unit_conversion = 1. #for now, keeping it 1.
	
	# ---- Output variables
	data_err = None
	list_headers = [None]
	
	# ---- Place holder for the output array
	Ncols = data[1].data.shape[1]
	data_err = np.zeros((len(data), 2, 2*pos_width+1, Ncols))
	data_err[:, 1, :, :] = infinity  # All data/errors are initially set to zero and infinity.
	
	for idx in range(1, len(data)):
	#for idx in range(98, 99):
		#Find the centroid along yaxis
		obpos = 'OBPOS' + str(idx)
		obposval = extr.header[obpos]
		obposval = int(obposval) #need to convert to int for slicing
		#print(f"idx: {idx}")
		# ---- Import data and unit conversion
		data_tmp = data[idx].data * unit_conversion
		err_tmp = err[idx].data * unit_conversion
		#print(data_tmp.shape)
		#print(data_tmp.shape)
		# ---- Trim the data
		#print(f"obposval: {obposval}")
		data_tmp = data_tmp[obposval-pos_width:obposval+pos_width+1]
		err_tmp = err_tmp[obposval-pos_width:obposval+pos_width+1]
		#print(data_tmp.shape)
		
		# ---- Apply preprocessing
		ibool = np.logical_or(np.isnan(err_tmp), np.isnan(data_tmp), err_tmp <=0.)
		data_tmp[ibool] = 0.
		err_tmp[ibool] = infinity
		
		# ---- Save data
		#print(data_err[idx, 0].shape)
		#print("------")
		data_err[idx, 0] = data_tmp
		data_err[idx, 1] = err_tmp
	
		# ---- Save header
		header_tmp = data[idx].header
		list_headers.append(header_tmp)
		
		#To debug
		#spectra_viewer(data_tmp)
		
	return data_err, list_headers

def bit_from_header(header):
	name = header["SLITOBJ"]
	if name == "stars":
		name = 2**1
	elif name == "gal":
		name = 2**2
	return int(name)

def extract_single_data(data_err, list_headers, specnum):
	"""
	Extract single spectrum data, err, header from the list, arr provided by Preprocessor.
	- specnum: Target object number in a file. Ranges from 1 through approx. 140.
	"""
	header = list_headers[specnum]
	data = data_err[specnum, 0]
	err = data_err[specnum, 1]
	
	return data, err, header

def ivar_from_err(err):
	return 1./np.square(err)

def naive_profile(data, ivar):
	"""
	Returns naive normalized profile
	"""
	K = np.sum(data * ivar, axis = 1) / np.sum(ivar, axis = 1)
	K /= np.sum(K) # Normalization
	return K

def produce_spec1D(data_err, list_headers, sig_K):
	"""
	Given 2D spectrum and the extraction kernel width sig_K,
	produce 1D spectra (Ntargets+1, 2, Ncols) and their inverse variance.
	"""
	data_ivar_1D = np.zeros((data_err.shape[0], 2, data_err.shape[3]))
	for specnum in range(1, len(list_headers)):
		data, err, header = extract_single_data(data_err, list_headers, specnum)
		ivar = ivar_from_err(err)

		spec1D_ivar = np.sum(np.square(K_T) * ivar, axis=0)
		spec1D = np.sum(K_T * data * ivar, axis=0) / spec1D_ivar
		
		data_ivar_1D[specnum, 0] = spec1D
		data_ivar_1D[specnum, 1] = spec1D_ivar
	return data_ivar_1D
		
def SIDE_from_header(header):
	return header["SIDE"]

def index_edges(data, num_thres=20):
	"""
	Given long postage stamp of data, return the edges.
	"""
	idx_min = 0
	idx_max = data.shape[1]-1
	tally = np.sum(data == 0., axis=0)
	while tally[idx_min] > num_thres:
		idx_min += 1
	while tally[idx_max] > num_thres:
		idx_max -=1
	return idx_min, idx_max

def gauss_fit2profile(K):
	# ---- Grid search for mu and sigma for the best Gaussian representation of the empirical kernel.
	Nrows = 32 
	mu_arr = np.arange(5, 20, 0.1)
	sig_arr = np.arange(1., 3., 0.05)
	chi_arr = np.zeros((mu_arr.size, sig_arr.size))
	x_arr = np.arange(0, Nrows, 1)
	for i, mu in enumerate(mu_arr):
		for j, sig in enumerate(sig_arr):
			A = np.exp(-np.square(x_arr-mu) / (2 * sig**2)) / (np.sqrt(2 * np.pi) * sig)
			chi_arr[i, j] = np.sum(np.square(K - A))
	# ---- Best fit
	idx = np.unravel_index(np.argmin(chi_arr), chi_arr.shape)
	mu_best = mu_arr[idx[0]] 
	sig_best = sig_arr[idx[1]]
	
	return mu_best, sig_best

def extract_stellar_profiles(masknumber, data_err, list_headers):
	K_collection = []
	
	count = 0 #counter for how many stellar profiles
	for specnum in range(1, len(list_headers)):
		data, err, header = extract_single_data(data_err, list_headers, specnum)
		ivar = ivar_from_err(err)

		BIT = bit_from_header(header)
		
		# ---- Perform optimal extraction 1D spectrum from 2D
		if (BIT == 2):
			K = naive_profile(data, ivar)
			K_collection.append(K) # Collect K into a list.
			
			#Plot stellar profile; optional for diagnosis
			crval1 = header['CRVAL1']
			cdelt1 = header['CDELT1']
			wavegrid = crval1 + cdelt1 * np.arange(data.shape[1])
			
			plt.imshow(data, extent = (crval1, wavegrid[-1], data.shape[0], 0), \
			cmap = 'gray', aspect = 'auto', vmin = 0, vmax = 1)
			plt.title("mask: " + str(masknumber) + " idx: " + str(specnum) +" count: " + str(count))
			
			#check if directory to save plots exist; if not, create one
			if not os.path.exists("../results/stellar_2D/" + str(masknumber) + "/"):
				os.makedirs("../results/stellar_2D/" + str(masknumber) + "/")
				
			plt.savefig("../results/stellar_2D/" + str(masknumber) + "/" \
			"idx-" + str(specnum) +"-count-" + str(count) + ".png",\
			dpi = 250, bbox_inches = 'tight')
			print("mask: " + str(masknumber) + " idx: " + str(specnum) +" count: " + str(count))
			count = count + 1
	return K_collection
	
def K_exp_profile(mu, sig, beta = 2, Nrows = 21): 
	"""
	Generate gaussian extraction profile of length Nrows
	given mu and sig.
	"""
	x_arr = np.arange(0, Nrows, 1)
	K_gauss = np.exp(-0.5 * (np.abs(x_arr - mu)/(sig))**beta)
	K_gauss /= np.sum(K_gauss)
	
	return K_gauss
		
def remove_outlier(arr, std_thres = 2):
	"""
	Remove outliers in 1D array by sigma clipping.
	"""
	std = np.std(arr)
	mu = np.median(arr)
	
	return arr[(arr - mu) < (std_thres * std)]


def extraction_kernel_sig(K_collection):
	"""
	Based on the extracted stellar profiles, 
	compute a reasonable gaussian extraction kernal
	width (sigma).
	"""
	# Array of estimates gaussian means and widths
	K_gauss_mus = np.zeros(len(K_collection))
	K_gauss_sig = np.zeros(len(K_collection))
	for i in range(len(K_collection)):
		mu_best, sig_best = gauss_fit2profile(K_collection[i])    
		K_gauss_mus[i] = mu_best
		K_gauss_sig[i] = sig_best

	return np.median(K_gauss_sig)

def K_gauss_profile(mu, sig, Nrows = 32):
	"""
	Generate gaussian extraction profile of length Nrows
	given mu and sig.
	"""
	x_arr = np.arange(0, Nrows, 1)
	K_gauss = np.exp(-(x_arr - mu)**2 / (2 * sig**2))
	K_gauss /= np.sum(K_gauss)
	
	return K_gauss

def plot_kernels(masknumber, K_collection, K_extract, fname):
	"""
	Plot the collection of stellar kernels and the ultimate
	extraction kernel at the center.
	Modify to include masknumber to plot pathological cases
	"""
	
	K_collection = np.array(K_collection)
	fig, ax = plt.subplots(1, figsize=(10, 5))
	for i in range(len(K_collection)):
		if(masknumber == 1624 or masknumber == 1625):
			if(K_collection[i, 0] == np.min(K_collection[:,0])):
				print(i)
				ax.plot(K_collection[i], c="green", lw=0.5)    
			else:
				ax.plot(K_collection[i], c="red", lw=0.5)
		elif(masknumber == 1627):
			if(K_collection[i, 0] == np.max(K_collection[:,0])):
				print(i)
				ax.plot(K_collection[i], c="green", lw=0.5)    
			else:
				ax.plot(K_collection[i], c="red", lw=0.5)
		else:
			ax.plot(K_collection[i], c="red", lw=0.5)
	ax.plot(K_extract, c="blue", lw=1.5)
	ax.axhline(y=0, c="black", ls="--", lw=1.)
	plt.savefig(fname, dpi=250, bbox_inches="tight")
	plt.close()

	return

def K_gauss_profile(mu, sig, Nrows = 32):
	"""
	Generate gaussian extraction profile of length Nrows
	given mu and sig.
	"""
	x_arr = np.arange(0, Nrows, 1)
	K_gauss = np.exp(-(x_arr - mu)**2 / (2 * sig**2))
	K_gauss /= np.sum(K_gauss)
	
	return K_gauss

def produce_spec1D(data_err, list_headers, K_extract, mu, fname_prefix=None, verbose=True):
	"""
	Given 2D spectrum and the extraction kernel K_extract,
	produce 1D spectra (Ntargets+1, 2, Ncols) and their inverse variance.
	"""
	data_ivar_1D = np.zeros((data_err.shape[0], 2, data_err.shape[3]))
	K_extract = K_extract[:,np.newaxis]
	for specnum in range(1, len(list_headers)):
		if verbose and ((specnum % 10) == 0): #Show progress every 10 spectra
			print("Processing spec num: %d" % specnum)
		data, err, header = extract_single_data(data_err, list_headers, specnum)
		ivar = ivar_from_err(err)

		spec1D_ivar = np.sum(np.square(K_extract) * ivar, axis=0)
		spec1D = np.sum(K_extract * data * ivar, axis=0) / spec1D_ivar
		
		#treat NaNs in spec1D data same as preprocess_bino and set them to 0.
		ibool = np.isnan(spec1D)
		spec1D[ibool] = 0.
		
		data_ivar_1D[specnum, 0] = spec1D
		data_ivar_1D[specnum, 1] = spec1D_ivar
		
		"""#extract edges along wavelength
		idx_min, idx_max = index_edges(data, num_thres=14)
		
		#produce wavegrid to see in spectra in range
		crval1, wavegrid = wave_grid(data, header)
			
		if fname_prefix is not None:
			plt.close()
			# ---- Spec figures
			
			fname = fname_prefix + "spec%d-2D.png" %specnum
			fig, ax = plt.subplots(1, figsize=(17, 1))
			ax.imshow(data, aspect="auto", cmap="gray", interpolation="none",\
			extent =  (crval1, wavegrid[-1], data.shape[0], 0), vmin=0., vmax=1)
			ax.set_xlim(wavegrid[idx_min], wavegrid[idx_max])
			ax.axhline(y=mu+0.5, c="red", ls="--", lw=0.8)
			ax.set_xlabel(r"$\lambda$ [nm]", fontsize = 15)
			plt.savefig(fname, dpi=250, bbox_inches="tight")
			plt.close()
			
			# ---- Histogram of centers determined
			fname = fname_prefix + "spec%d-centers.png" %specnum
			fig, ax = plt.subplots(1, figsize=(7, 3))
			ax.hist(row_centers, bins=np.arange(0.5, 32.5, 1), histtype="step", color="black", normed=True)
			ax.plot(K_extract, c="red", label="K_stellar")
			ax.axvline(x=mu, c="red", ls="--", lw=0.4)
			plt.savefig(fname, dpi=200, bbox_inches="tight")
			plt.close()"""

	return data_ivar_1D

def wavegrid_from_header(header, Ncols):
    """
    Construct a linear grid based on the header
    and a user specified number of columns.
    """
    x0 = header["CRVAL1"] * 10
    dx = header["CDELT1"] * 10
    return x0 + np.arange(0, Ncols, 1.) * dx
    
def find_nearest_idx(arr, x):
	return np.argmin(np.abs(arr-x))