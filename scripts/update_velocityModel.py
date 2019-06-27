"""This script looks at excel notebook, recalculates the 
velocity and amplitude based on the correct redshift.

This requires a modification of SNR_calculator() which ONLY
searches over the width space, and NOT the redshift space."""

from utils import datareader, wave_grid, Window, Model
from global_var import *
import numpy as np
import os
import matplotlib.pyplot as plt
from time import time

#user input of maskname; as str
maskname = sys.argv[1]
flag = (sys.argv[2].lower() == 'true') #user entered input whether bluer or redder mask

if(flag == True):
	maskname_b = maskname
	maskname_r = str(int(maskname) + 1)
else:
	maskname_r = maskname
	maskname_b = str(int(maskname) - 1)

data = datareader(maskname)

def SNR_calculator_modified(maskname, data, z, rel_strngth, fudge_factor):
	"""modified version of SNR_calculator"""
	
	#Read data
	image = data['data_ivar'][:, 0, :]
	ivar = data['data_ivar'][:, 1, :]
	crval1, wg = wave_grid(data)
	
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
	SNRs = np.zeros((image.shape[0], sigma_v.size))
	
	#Save all the amplitudes and chi_sq to pass this to the PeakZoom function
	#Amps -> (z x num. of objects x width)
	#del_chi_sq -> (z x num. of objects x width)
	Amps = np.zeros(SNRs.shape) 
	del_chi_sq = np.zeros(SNRs.shape)
	
	#Store all the widths b/c widths are a func. of z
	widths = np.zeros((sigma_v.size))
	
	#-------------------------------------------#
	#WORK ON THIS LATER
	#Save all the medians to pass to PeakZoom function
	medians = np.zeros((image.shape[0]))
	#-------------------------------------------#
		
	wg2 = Window(z, wg)
	widths = widthlist(z) #Annoying bug. MUST save widths becuase widths = widths(z). Previously was using last z 
	model = Model(z, wg2, widths, relative_strength = rel_strngth)
		
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
	medians = median_val
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
	
	SNRs = SNR
	Amps = Amp
		
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
	del_chi_sq = del_chi2
	
	SNRs_final = SNRs #This maintains the indices
	Amps_final = Amps
	del_chi_sq_final = del_chi_sq
	
	return widths, SNRs_final, Amps_final, del_chi_sq_final, medians
	

start = time()
widths, SNRdata, Ampdata, delChi2data, mediansdata = \
SNR_calculator_modified(maskname, data, z, relative_strength, fudge_factor)
end = time()

print(f"total time: {end-start}")