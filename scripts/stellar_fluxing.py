"""This script produces stellar fluxing plots given a mask"""

import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
import json
from global_var import *
import sys
import os
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from scipy.signal import medfilt
from utils import datareader, wave_grid
from sedpy.observate import load_filters, getSED
from calbino import choose_model, get_library

plt.style.use("fivethirtyeight")

#user input of maskname; as str
maskname = sys.argv[1]

#Read in spec1d npz datafile; data['headers'] and data['data_ivar']
data = datareader(masknumber)

# find all the standard stars
id_array = []
ra_array = []
dec_array = []
for i in range(1,len(data['headers'])):
    if((data['headers'][i]['SLITOBJ'] == 'stars') | (data['headers'][i]['SLITOBJ'] == '2')):
        print(i)
        id_array.append(i)
        ra_array.append(data['headers'][i]['SLITRA'])
        dec_array.append(data['headers'][i]['SLITDEC'])
        
id_array = np.asarray(id_array)
ra_array = np.asarray(ra_array)
dec_array = np.asarray(dec_array)

#tabulate all the standard stars in the mask
std_stars_df = pd.DataFrame(id_array.astype(np.int), columns=['Slit ID'])
std_stars_df['RA'] = ra_array
std_stars_df['Dec'] = dec_array

##---Cross match catalogue with SDSS for colour---##
deg2arcsec=3600

def match_cat1_to_cat2(ra1, dec1, ra2, dec2):
    """
    "c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)  
    catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)  

    idx are indices into catalog that are the closest objects to
    each of the coordinates in c, 
    d2d are the on-sky distances between them, and d3d are the
    3-dimensional distances."
    -- astropy documentation.  

    Fore more information: 
    http://docs.astropy.org/en/stable/coordinates/matchsep.html#astropy-coordinates-matching 
    """    
    cat1 = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)  
    cat2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  
    idx, d2d, d3d = cat1.match_to_catalog_sky(cat2)
    
    return idx, d2d.degree

def crossmatch_cat1_to_cat2(ra1, dec1, ra2, dec2, tol=1./(deg2arcsec+1e-12)):
    """
    Return indices of cat1 (e.g., DR3) and cat2 (e.g., DEE2) cross matched to tolerance. 

    Note: Function used to cross-match DEEP2 and DR3 catalogs in each field 
    and test for any astrometric discrepancies. That is, for every object in 
    DR3, find the nearest object in DEEP2. For each DEEP2 object matched, 
    pick DR3 object that is the closest. The surviving objects after these 
    matching process are the cross-matched set.
    """
    
    # Match cat1 to cat2 using astropy functions.
    idx_cat1_to_cat2, d2d = match_cat1_to_cat2(ra1, dec1, ra2, dec2)
    
    # Indicies of unique cat2 objects that were matched.
    cat2matched = np.unique(idx_cat1_to_cat2)
    
    # For each cat2 object matched, pick cat1 object that is the closest. 
    # Skip if the closest objects more than tol distance away.
    idx1 = [] # Place holder for indices
    idx2 = []
    tag = np.arange(ra1.size,dtype=int)
    for e in cat2matched:
        ibool = (idx_cat1_to_cat2==e)
        candidates = tag[ibool]
        dist2candidates = d2d[ibool]
        # Index of the minimum distance cat1 object
        if dist2candidates.min()<tol:
            idx1.append(candidates[np.argmin(dist2candidates)])
            idx2.append(e)
    
    # Turning list of indices into numpy arrays.
    idx1 = np.asarray(idx1)
    idx2 = np.asarray(idx2)
    
    # Return the indices of cat1 and cat2 of cross-matched objects.
    return idx1, idx2
	
#need to convert from pandas to astropy Table for astropy operations
std_stars_df = Table.from_pandas(std_stars_df)
print("Found number of standard stars: ", len(std_stars_df))

#read in standard bino file
std_file = pd.read_csv("../data/binostandardfiles/standard_bino.r0_45.dm6_p6.txt", \
                      sep = " ")
std_file = Table.from_pandas(std_file)

#first indices correspond to first catalog, second with second catalog
idxcat1, idxcat2 = crossmatch_cat1_to_cat2(std_stars_df['RA'], std_stars_df['Dec'], \
                        std_file['ra'], std_file['dec'])
						
#sift out stars in mask
mask_std_stars = std_file[idxcat2]
std_stars_df = std_stars_df[idxcat1]

len("Found number of standard stars in the SDSS file: ", len(mask_std_stars))

##---Fluxing---##
#in this section, generate spectra of desired standard stars
image = data['data_ivar'][:, 0, :]
crval1, wg1 = wave_grid(data)

def spec_fetcher(idx, trim = 100):
    """fetches spectra of desired standard star
    Parameters
    ----------
    idx: index of the star in the given mask
    trim: amount of pixels trimmed from the end of the spectra"""
    
    imagetmp = image[idx, :-trim]
    wg = wg1[:-trim]
    
    #plot spectra
    #plt.plot(wg, imagetmp)
    
    return wg, imagetmp

def standard_star_compare():
	# This is the filters and magnitudes for the calibrator star
	filters = load_filters(["sdss_g0", "sdss_r0", "sdss_i0", "sdss_z0"])
	
	# Get a reasonable set of model spectra
	libwave, libflux, libparams = get_library()
	
	
	#NEED TO CHANGE 0 TO INDEX BASED ON CROSSMATCH
	#loop over 0
	for i, val in enumerate(std_stars_df):
		
		star_mags = np.array([mask_std_stars['psfMag_g'][i],\
						mask_std_stars['psfMag_r'][i],\
						mask_std_stars['psfMag_i'][i],\
						mask_std_stars['psfMag_z'][i]])
		
		#bino spectrum for calibration star
		data_wave, data_flux = spec_fetcher(val['Slit ID'])
		data_flux *= 1e-19  # header units
		
		# choose the model with colors closest to the calibrator
		best_model = choose_model(star_mags, filters, libwave, libflux)
		
		# Now work out the normalization of the model from the average magnitude offset
		best_sed = getSED(libwave, best_model, filters)
		dm = np.mean(star_mags - best_sed)
		conv = 10**(-0.4 * dm)
		
		# Here, finally, is the fluxed model (erg/s/cm^2/AA)
		fluxed_model = best_model * conv
		
		# Now get the model on the same wavelength vector as the data
		z = 0.0 # redshift of the star, if known. #?
		a = (1 + z)
		fluxed_model_interp = np.interp(data_wave, libwave * a, fluxed_model)
		calibration = data_flux / fluxed_model_interp
		
		# You probably want to median filter the calibration vector.  # Perhaps
		# after some sigma clipping.  differences on small scales could be due to
		# model imperfections (wrong metallicity, wrong gravity for model, LSF
		# mismatch)
		# you could also fit the calibration vector with a polynomial, taking into
		# account errors
		smoothed_calibration = medfilt(calibration, 101)
		
		#fig, axes = plt.subplots(1, 1, sharex=True, figsize=(13, 11))
		#plt.plot(data_wave, calibration, label="raw calibration")
		plt.plot(data_wave, smoothed_calibration, label="smoothed calibration")
		plt.xlabel("Wavelength [A]")
		plt.ylabel("actual / input")
		plt.savefig("../results/stellar_fluxing/" + maskname + ".png", \
		dpi = 250, bbox_inches = 'tight')
		#print(np.min(smoothed_calibration))
		#ax.legend()
		#ax.set_ylabel("actual / input")
		#plt.ylim([0, 0.1])
		#plt.xlim([7500, 7700])
		
standard_star_compare()