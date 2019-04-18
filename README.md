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
of slits for a given mask `python produce_spec1d.py maskname`.
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
other products that will be used to measure redshifts `python SNR_Amp.py maskname`.
	* It will require a user input of the masknumber to generate all the outputs.
	* It will generate the following data products, where * stands for the number 
	denoting [O II] doublet ratio of the 1st line (3727) wrt the 2nd line (3729):
		* /results/outputdata/masknumber/masknumber-*-Amp.npy  -- **Amplitude** of [O II] 
		doublet hypotheses (slit index, redz, width)
		* /results/outputdata/masknumber/masknumber-*-delChi2.npy  -- **delta Chi2** 
		of existence of [O II] doublet vs no doublet (slit index, redz, width)
		* /results/outputdata/masknumber/masknumber-*-SNR.npy  -- **SNR** of [O II] doublet 
		hypotheses (slit index, redz, width)
		* /results/outputdata/masknumber/masknumber-*-widths.npy  -- **widths** of [O II] 
		doublet hypotheses (slit index, redz)
		* /results/outputdata/masknumber/masknumber-*-widths.npy   -- **z_range** list of 
		redshift hypotheses (slit index)

3. Next, run /scripts/**best_model.py** to generate a list of best models, i.e. redshifts 
and widths, that best explain the data for a given mask. If a mask has *n* slits, then there
will be *n* redz and widths, each corresponding to exactly one slit.
	* This script takes the outputs of SNR_Amp.py as input. The user provides the maskname 
	for which the models will be generated.
	* It will generate the following data products:
		* /results/Max_z_n_width/maskname.txt which is a csv file of three columns: 
		slit index, redz and width
		* /results/Max_z_n_width/maskname-indices.npy which is a npy file containing
		the idx of best redz and width for every slit for a given mask. The idx here
		corresponds to the idx for *SNR.npy , *delChi2.npy and *Amp.npy produced by 
		SNR_Amp.py. We produce this as it will be later used by two additional 
		/scripts/scripts peak_zoom.py and /scripts/hyp_test.py.
		* /results/histograms/redshift/maskname.png which is the initial redshift 
		histogram.
		* /results/histograms/width/maskname.png which is the initial velocity dispersion 
		histogram.
		
4. Next, run /scripts/**hyp_test.py** to generate hypothesis plots which tries to 
find all the other associated lines for a given emission lines, i.e. if we assume that 
the line we are seeing is a [O II] doublet, the produced plot will find the postage 
stamps where we should expect to see Hbeta, [O III] 4960, [O III] 5007 and Halpha.
	* This script takes user provided maskname and a "True" or False (or anything else) 
	binary as a boolean. The boolean denotes whether the mask in question is the bluer mask
	or the redder one. Since each region of the sky was observed in the blue wavelength
	range (~4000 - 8200 A) and red wavelength range (~7000 - 10500 A), we use the boolean 
	to determine whether the mask in question is bluer or redder. By knowing the boolean,
	the script can then understand where it should search for additional emission lines.
	* It will generate the following data products:
		* /results/hyp_test/maskname/maskname-Slit-*.png where * denotes the slit index number. 
		There will be a total of len(mask) number of plots in every folder.
		
5. Next, run /scripts/**excel_merger.py** to construct an Excel file in /results/Max_z_n_width
 that will combine outputs of *best_model.py* for both the blue and the red mask into an 
 excel file that can be used to document any false positives. The format of the file name 
 is *bluemask+redmask_visual-inspected.xlsx*.
	* The constructed file is then used in conjunction with the data product of *hyp_test.py*
	to check for any false positives. If a false positive is detected, then the correct 
	redshift and width is noted in the Excel file.
	* **If the false positive is a clear artifact, then the associated redshift and width 
	values should be deleted from the excel file.**
	* In addition to identifying false positives, this step also assigns confidence to
	redshift measurement. A thorough description of what each confidence class looks like 
	is described in /results/Max_z_n_width/**visual-inspection-keywords.readme**.
	* There are two kinds of false positives - *misclassified line* and *artifacts appearing 
	as signal*. Each of these false positives are treated differently. Point 6 and 7 
	deals with these two cases.
	
6. When inspecting the outputs of *hyp_test.py*, **always first look for misclassified lines**. 
This is because the script **redz_line_calculator.py** can update the excel sheet that is 
output by **excel_merger.py** with the proper redshift, confidence level and notes. 
**Remember that you should do this step while keeping the excel sheet closed otherwise 
redz_line_calculator.py cannot access the file**. A description of what Notes need to be 
provided in this step is given in /results/Max_z_n_width/visual-inspection-keywords.readme.
	* The input parameters of *redz_line_calculator.py* is somewhat complicated. 
	In order to update the excel sheet, enter the following parameters in the command line: 
	`python redz_line_calculator.py bluemask redmask flag confidence ion idx`
	* bluemask and redmask represent the four digit numerical codes of the mask as
	found in docs/folder_to_mask.txt.
	* flag refers to whether the column we want to update is the blue mask or the 
	red mask one. True if blue mask, other if red.
	* confidence denotes the confidence level as described in step 5.
	* ion refers to what the user suspects the falsely identified [O II] line to actually be.
	For example, if the falsely identified line is identified to be [O III] 5007, then the user
	inputs o2. In total, there are four lines that the user can propose for which the redshifts
	can be updated -- h beta, [O III] 4960, [O III] 5007 and h alpha.
	* idx refers to the slit index.
	
7. After checking for misclassified lines, check for artifacts. These are falsely identified 
[O II] lines that are due to cosmic rays, bad pixels, edge effect, etc. To do this, 
**run hyp_additional.py**.
	* It requires the following input `python hyp_additional.py maskname flag slitidx nthidx`:
		* maskname
		* flag -- whether mask is blue mask or red mask
		* slitidx
		* nthidx -- if the idx corresponding to the lowest deltaChi2 is not the true [O II]
		line, then out a hypothesis plot similar to the output of *hyp_test.py* where nthidx
		denotes the nth lowest deltaChi2. 
		* Notice that this script requires some development as it is part science part art.
		The code can tell you whether your proposed nthidx hypothesis redshift is closer than
		0.01 to the best hypothesis and whether the SNR for nthidx is less than SNR = 10, in 
		which case we ignore those answers. So sometimes you will have to do trial-and-error 
		to find actual lines. 












