""""This script generates SNR and Amp datafiles for the masks. The outputs
can then be used to generate visualizations."""

from utils import datareader, SNR_calculator
import numpy as np
import sys
import os 

maskname = sys.argv[1]
	
data = datareader(maskname)

relative_strength = 0.7
fudge_factor = 1.

#velocity dispersion array
sigma_v = np.arange(0, 201, 10)

#Time the code
from time import time
start = time()

#Calculate SNR, Amp, delChi2
z_range, widths, SNRdata, Ampdata, delChi2data, mediansdata = \
SNR_calculator(maskname, data, rel_strngth = relative_strength, fudge_factor = fudge_factor)

end = time()

print(f"total time: {end-start}")

#check if directory to save plots exist; if not, create one
if not os.path.exists("../results/outputdata/" + maskname + "/"):
	os.makedirs("../results/outputdata/" + maskname + "/")

np.save("../results/outputdata/" + maskname + "/" + maskname + "-" + str(relative_strength) + "-SNR.npy", SNRdata)
np.save("../results/outputdata/" + maskname + "/" + maskname + "-" + str(relative_strength) + "-Amp.npy", Ampdata)
np.save("../results/outputdata/" + maskname + "/" + maskname + "-" + str(relative_strength) + "-delChi2.npy", delChi2data)
np.save("../results/outputdata/" + maskname + "/" + maskname + "-" + str(relative_strength) + "-z_range.npy", z_range)
np.save("../results/outputdata/" + maskname + "/" + maskname + "-" + str(relative_strength) + "-widths.npy", widths)