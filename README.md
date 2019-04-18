# DESI Emission-Line Galaxies Target Selection Validation
### Author: Tanveer Karim

The goal of this project is to measure the efficiency of various target
selection strategies for emission-line galaxies (ELGs) within the Dark
Energy Spectroscopic Instrument (DESI) experiment. DESI is a multi-year
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

## Running the Pipeline
1. Download data from the MMTO server. The folder should have the following files:
   *obj_abs_err_slits_lin.fits
   *obj_abs_slits_extr.fits
   *obj_abs_slits_lin.fits
   *obj_counts_qlook.fits
	
2. 
In Binospec_DataProcessor/Tanveer/Specz_Extractor/scripts, there should exist a file called utils_spec1d.py that contains all the necessary functions to run this pipeline. Update this script if you want to update the pipeline. In addition, there is a file called global_var.py that contains the global variables used throughout all the scripts. Import this file too to run the scripts.
