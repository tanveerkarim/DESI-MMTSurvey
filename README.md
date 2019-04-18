# DESI Emission-Line Galaxies Target Selection Validation
### Author: [Tanveer Karim](tanveerkarim.com)

The goal of this project is to measure the efficiency of various target
selection strategies for emission-line galaxies (ELGs) within the [Dark
Energy Spectroscopic Instrument (DESI)](https://www.desi.lbl.gov/) experiment. DESI is a multi-year
dark energy experiment whose goal is to study the expansion of the 
Universe over the past 12 billion years and create the most precise 
3D map of the Universe to-date. In order to do this, DESI will be observing
30 million galaxies, 17 million of which are ELGs. 

This project is the data analysis pipeline for the target selection validation 
survey. The survey was conducted using the [MMT Binospec instrument](https://www.cfa.harvard.edu/mmti/binospec.html).
15 high-target regions were chosen for this validation survey and the various
target selection strategies proposed a total of ~3,000 galaxies to be observed.
With the help of this pipeline, we check whether these galaxies meet the redshift
and [O II] flux threshold of DESI, and which strategies performed the best.

## Getting Started
### Prerequisites

### Assumptions about Data
1. Each of the regions have a blue mask and a red mask.
2. The blue mask and red mask data are in separate folders with all the requisite files.
3. The blue masks and red masks are sequential, i.e. redmask = bluemask +1

### Data Storage and Important Functions
1. All the data should be stored in data/data_folders/. Download data from the MMTO server. 
Keep in mind the assumptions explained above. The folder should have the following files:
   * obj_abs_err_slits_lin.fits
   * obj_abs_slits_extr.fits
   * obj_abs_slits_lin.fits
   * obj_counts_qlook.fits
 
2. /scripts/utils.py, contains all the necessary functions to run this pipeline. 
The original version was written by Jae Lee and the subsequent versions have been 
maintained by the current author. Update this script if you want to update the 
pipeline. In addition, /scripts/global_var.py contains the global variables used 
throughout all the scripts. Import this file too to run the scripts.

## Running the Pipeline
1. *First*, run /scripts/**produce_spec1d.py** to produce 1d spectra 
of slits for a given mask.
	* It will require a user input of the masknumber to generate all the outputs. 
	masknumber is a four digit code for a mask produced by Binospec. The dictionary
	that contains the information of masknumber to actual maskname is in 
	docs/folder_to_mask.txt.
	* Set the pos_width to 7 to extract a 15 pixel wide window where the 
	object will be in the middle of the window.
	* Make sure that the directory of data_dir in utils.py is pointing 
	to the data directory.
	* This code will generate terminal outputs that tells the user which mask 
	is being analysed, what idx correspond to standard stars and will produce
	a count for how many stars are there. It will also show the overall processing 
	of all the spectra in steps of 10.
	* In addition, this code will generate three outputs that are saved in 
	various directories:
		* /results/kernels/masknumber-kernel.png shows all the standard star 
		kernels and the average kernel that will be used for extraction. Note 
		that rather than using a Gaussian, we are using a ![equation](https://latex.codecogs.com/gif.latex?%5Cexp%5Cleft%28%20-0.5%5Cleft%28%20%5Cfrac%7Bx%20-%20%5Cmu%7D%7B%5Csigma%7D%20%5Cright%29%5E4%20%5Cright%29)
		function; hence the average extraction kernel is fatter.
		* /results/stellar_2D/maskname/idx-*-count-*.png contains all the 2D 
		spectra of the standard stars where idx number refers to slit number and 
		count number refers to the ordinal number given to the standard star as 
		explained in the previous point.
		* data/npz_files/*-spec1d.npz contains masknumber-spec1d.npz files that contains 
		all the 1D spectra based on the extraction kernel and 1D headers in the 
		form of a npz file.
		
2. Next, run /scripts/**SNR_Amp.py** to generate the matrix of signal-to-noise and 
other products that will be used to measure redshifts.
	* It will require a user input of the masknumber to generate all the outputs.
	* It will generate the following data products, where * stands for the number 
	denoting [O II] doublet ratio of the 1st line (3727) wrt the 2nd line (3729):
		* /results/outputdata/masknumber/masknumber-*-Amp.npy  -- **Amplitude** of [O II] doublet hypotheses (slit index, redz, width)
		* /results/outputdata/masknumber/masknumber-*-delChi2.npy  -- **delta Chi2** of existence of [O II] doublet vs no doublet (slit index, redz, width)
		* /results/outputdata/masknumber/masknumber-*-SNR.npy  -- **SNR** of [O II] doublet hypotheses (slit index, redz, width)
		* /results/outputdata/masknumber/masknumber-*-widths.npy  -- **widths** of [O II] doublet hypotheses (slit index, redz)
		* /results/outputdata/masknumber/masknumber-*-widths.npy   -- **z_range** list of redshift hypotheses (slit index)





