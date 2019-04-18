import numpy as np

#preprocess_bino
pos_width_1626 = 7 #1626 is smaller in width so need 7; others work with 10
pos_width = 7

#[O II] doublet
lambda27 = 3727.092
lambda29 = 3729.875
separation = (lambda29 - lambda27)/2 #separation between the emission lines
lambda0 = lambda27 + separation #Midpoint of the gaussian emission lines in restframe

#velocity dispersion array
sigma_v = np.arange(0, 201, 5)

#Other lines; values from SDSS
OII = lambda0
hB0 = 4862.68
OIII1 = 4960.295
OIII2 = 5008.240
ha = 6564.61   

#plotting variables
alpha = 0.7
vmin = -0.5; vmax = 1 #for imshow
fontsize = 25
font = 'serif'
aspect = 'auto'

#functional parameters for [O II] modelling
relative_strength = 0.7 #relative strength ratio between left line and right line
fudge_factor = 1. #multiplicative factor over slit width

#for mask_summary.py histogram plot
histtype = 'step'
lw = 3 #line width of histogram
#range of histogram for same binning
hist_left = 0.
hist_right = 1.6