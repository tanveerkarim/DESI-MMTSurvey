"""This script generates stellar kernels and produces 1D spectra for a given mask"""

from utils import *
from global_var import *
import sys

#user input of maskname; as str
#User inputs in sequential order: maskname
maskname = sys.argv[1]
masknumber = int(maskname)

#pos_width = 7 #based on HSC_Analysis/kernels/ 15 pixels total seems to encapsulate all the values
data_err, list_headers = preprocess_bino(data_dir = "../../../../DATA_MAY18/" + maskname + "/", \
	masknumber = masknumber, pos_width = pos_width)

# ---- Compute extraction kernel width
K_collection = extract_stellar_profiles(masknumber, data_err, list_headers)
plate_scale_x = .25 #arcsec/pixel
size = 1.5**2*np.pi/4 #See Daniel notes
sig_extract = (size/np.sqrt(2*np.pi))/plate_scale_x #default sigma = size/sqrt(2*pi) assuming 1.7 tophat"; converted to px
	
mu = pos_width + 0.5
	
K_extract = K_exp_profile(mu, sig_extract, beta = 4, Nrows = 2*pos_width+1)
plot_kernels(masknumber, K_collection, K_extract, "../results/kernels/"+str(masknumber)+"-kernel.png")

# ---- Perform 1D extraction
data_ivar_1D = produce_spec1D(data_err, list_headers, K_extract, mu)#, fname_prefix=fname_prefix)
np.savez("../npz_files/" + str(masknumber) +"-spec1d.npz", data_ivar = data_ivar_1D, headers = list_headers)
print("---------")

