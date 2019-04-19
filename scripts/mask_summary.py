"""This script generates redz summary data based on algorithms"""

import numpy as np
import pandas as pd
from utils import datareader
from global_var import *
import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys
import os
plt.style.use("fivethirtyeight")

def bit_from_header(header):
	name = header
	if name == "stars":
		name = 2**1
	elif name == "gal":
		name = 2**2
	return int(name)

##--ALL THE DICTIONARIES USED IN THIS CODE--##

#Dictionary in which subdataframes by confidence level only are stored
confidence_df_dict = {}

#Dictionary consisting of algorithm name and bit values 
bit_dict = {'fdr': 2**3, 'ndm': 2**4 + 2**6 + 2**7, 
 'rf': 2**8 + 2**9 + 2**10}

#Dictionary consisting of 9 dataframes by confidence class
#and algorithm class
conf_alg_df_dict = {}

#Dictionary consisting of confidence level naming convention;
#the keys are designed to iterate over for plotting purposes
#such as the histogram
conf_name_dict = {1:'3', 2:'3+2', 3:'3+2+1'}

#Dictionary for storing all the statistics
statistics_dict = {}
##--END OF LIST OF DICTIONARIES--##

#user input of maskname; as str
maskname_b = sys.argv[1]
maskname_r = sys.argv[2]

#read in data based on maskname
data_b = datareader(maskname_b)
data_r = datareader(maskname_r)

#Read in visually inspected xlsx file
df = pd.read_excel("../results/Max_z_n_width/" + maskname_b + "+" + maskname_r + "_visual-inspected.xlsx", header=0)

#Check whether blue mask redz agree with red mask redz
df.plot('z_b', 'z_r', marker = 'o', linestyle = "")
x = np.arange(0, 1.4, .1)
plt.plot(x, x, alpha = alpha, linestyle = "--")
plt.xlabel(r"$z_b$", fontsize = 15)
plt.ylabel(r"$z_r$", fontsize = 15)
plt.legend([])
plt.title(maskname_b + " vs " + maskname_r + " redshift comparison", family = 'serif', fontsize = 15)
#check if summary_statistics/ folder exists
output_dir = '../results/summary_statistics/'
if not os.path.exists(output_dir):
	os.makedirs(output_dir)
figname = output_dir + maskname_b + '+' + maskname_r + '_redz_compare.png'
plt.savefig(figname, dpi = 250, bbox_inches = 'tight')
plt.close()


#Store bit values in the dataframes
tmpbit = []
for i in range(1,len(df)+1):
    bit = bit_from_header(data_b['headers'][i]['SLITOBJ'])
    tmpbit.append(bit)
df['bitval'] = tmpbit

#Create subdataframes based on confidence levels:
confidence_df_dict['df_conf3'] = df.loc[((df['Confidence_b'] == 3) | (df['Confidence_r'] == 3))] #Confidence 3
confidence_df_dict['df_conf32'] = df.loc[((df['Confidence_b'] == 3) | (df['Confidence_r'] == 3)) | \
 ((df['Confidence_b'] == 2) | (df['Confidence_r'] == 2))] #Confidence 3+2
confidence_df_dict['df_conf321'] = df.loc[((df['Confidence_b'] == 3) | (df['Confidence_r'] == 3)) | \
 ((df['Confidence_b'] == 2) | (df['Confidence_r'] == 2)) | \
 ((df['Confidence_b'] == 1) | (df['Confidence_r'] == 1))] #Confidence 3+2+1

#calculate mean z to get the final z that will be used for statistics
for key in confidence_df_dict:
    confidence_df_dict[key]['z_final'] = confidence_df_dict[key][['z_b', 'z_r']].mean(axis = 1)
	
#Split the 3 confidence levels into 3 further levels on the basis of algorithms
for bitkey in bit_dict:
    for dfkey in confidence_df_dict:
        conf_alg_df_dict[dfkey + '_' + bitkey] = \
        confidence_df_dict[dfkey][np.bitwise_and(confidence_df_dict[dfkey]['bitval'], \
                                                 bit_dict[bitkey]) > 0]

##--HISTOGRAM PLOTTING--##										 
#Regular histogram
fig = plt.figure(1, figsize=(20,7))
gridspec.GridSpec(1,3)

for i in range(1,4):
    ax = plt.subplot2grid((1,3), (0, i-1))
    confidence_df_dict['df_conf' + "".join(conf_name_dict[i].split("+"))].hist('z_final', bins = 10, \
                                                             ax = ax, alpha = alpha, \
                                                            label = "Total", range = (hist_left, hist_right))
    
    #"".join([str1, str2]) concatenates all the strings into str1str2
    for bitkey in bit_dict:
        dfname = 'df_conf' + "".join(conf_name_dict[i].split("+")) + '_' + bitkey
        
        if(bitkey == 'rf'):
            conf_alg_df_dict[dfname].hist('z_final', bins = 10, histtype = histtype, \
                                         linewidth = lw, ax = ax, \
                                          label = bitkey.upper(), color = 'k', range = (hist_left, hist_right))
        else:
            conf_alg_df_dict[dfname].hist('z_final', bins = 10, histtype = histtype, \
                                         linewidth = lw, ax = ax, \
                                          label = bitkey.upper(), range = (hist_left, hist_right))
        ax.set_title('Confidence ' + conf_name_dict[i] + ' Redshifts', \
                     family = font, fontsize = fontsize)
        ax.legend(loc = 'best')
		
figname = '../results/summary_statistics/' + maskname_b + '+' + maskname_r + '_hist.png'
plt.savefig(figname, dpi = 250, bbox_inches = 'tight')
plt.close()

#Normalized histogram
"""fig = plt.figure(1, figsize=(20,7))
gridspec.GridSpec(1,3)

for i in range(1,4):
    ax = plt.subplot2grid((1,3), (0, i-1))
    confidence_df_dict['df_conf' + "".join(conf_name_dict[i].split("+"))].hist('z_final', bins = 10, \
                                                             ax = ax, density = True, alpha = alpha, \
                                                            label = "Total")
    
    #"".join([str1, str2]) concatenates all the strings into str1str2
    for bitkey in bit_dict:
        dfname = 'df_conf' + "".join(conf_name_dict[i].split("+")) + '_' + bitkey
        
        if(bitkey == 'rf'):
            conf_alg_df_dict[dfname].hist('z_final', bins = 10, histtype = histtype, \
                                         linewidth = lw, ax = ax, density = True, \
                                          label = bitkey.upper(), color = 'k')
        else:
            conf_alg_df_dict[dfname].hist('z_final', bins = 10, histtype = histtype, \
                                         linewidth = lw, ax = ax, density = True, \
                                          label = bitkey.upper())
        ax.set_title('Confidence ' + conf_name_dict[i] + ' Redshifts', \
                     family = font, fontsize = fontsize)
        ax.legend(loc = 'best')
		
figname = '../results/summary_statistics/' + maskname_b + '+' + maskname_r + '_hist_normalized.png'
plt.savefig(figname, dpi = 250, bbox_inches = 'tight')
plt.close()"""
##--END OF HISTOGRAM PLOTTING--##

##--SUMMARY STATISTICS TABLE--##
#calculates summary statistics based on dataframe and algorithm
def summary_stats(dataframe):
    """Returns summary statistics based on given dataframe
    Parameters
    ----------
    dataframe: confidence+algorithm dataframe for which stats
    is required; from conf_alg_df_dict
    """
    
    #Percentage of objects proposed by this confidence+algo in this mask
    nobj = len(conf_alg_df_dict[dataframe])
    percentage_in_mask = nobj/len(df)*100 #df is the master table
    
    #total number of objects in this algorithm in this mask
    #dataframe.split("_")[-1] finds whether last chunk in dataframe
    #name is fdr, ndm or rf. Then bit_dict calls the bit value
    nobjtot = len(df[np.bitwise_and(df['bitval'], bit_dict[dataframe.split("_")[-1]]) > 0.])
    
    #total success rate in alg
    totalz_success_rate_in_alg = nobj / nobjtot*100
    
    #highz = [1.1 +]; medz = [0.6, 1.1]
    highz_success_rate_in_alg = \
    len(conf_alg_df_dict[dataframe].loc[conf_alg_df_dict[dataframe]['z_final'] > 1.1]) / \
                                                                        nobjtot*100
    medz_success_rate_in_alg = \
    len(conf_alg_df_dict[dataframe].loc[(conf_alg_df_dict[dataframe]['z_final'] > 0.6) \
                        & (conf_alg_df_dict[dataframe]['z_final'] < 1.1)])/nobjtot*100
    
    #round all rates to 2 decimal places
    percentage_in_mask = np.round(percentage_in_mask, 2)
    totalz_success_rate_in_alg = np.round(totalz_success_rate_in_alg, 2)
    highz_success_rate_in_alg = np.round(highz_success_rate_in_alg, 2)
    medz_success_rate_in_alg = np.round(medz_success_rate_in_alg, 2)
    
    return nobj, percentage_in_mask, \
     totalz_success_rate_in_alg, highz_success_rate_in_alg, medz_success_rate_in_alg

def noZ_stats(dataframe):
    """Returns summary statistics of noZ class
    based on given dataframe.
    Parameters
    ----------
    dataframe: confidence+algorithm dataframe for which stats
    is required; from conf_alg_df_dict
    """
    
    #Total no. of objects in mask
    nobjtot_in_mask = len(df)
    
    #Total no. of objects in mask for a given algorithm
    nobjtot_in_alg = len(df[np.bitwise_and(df['bitval'], bit_dict[dataframe.split("_")[-1]]) > 0])
    
    #Total no. of noZ objects in given algorithm
    noZobj_in_alg = nobjtot_in_alg - len(conf_alg_df_dict[dataframe]) #dataframe must be of type 321
    noZ_percentage_in_alg = noZobj_in_alg/nobjtot_in_alg*100 #percentage in alg
    noZ_percentage_in_mask = noZobj_in_alg/nobjtot_in_mask*100 #percentage in mask
    
    return noZobj_in_alg, \
     np.round(noZ_percentage_in_mask,2), np.round(noZ_percentage_in_alg,2), '-', '-'

#add 'NoZ' as the 4th conf_name_dict val
conf_name_dict[4] = 'noZ'

#Calculate all the necessary statistics and store in statistics_dict
for bitkey in bit_dict:
    for i in range(1,5):
        #print(conf_name_dict[i])
        if(conf_name_dict[i] == 'noZ'):
            key = 'df_conf321_' + bitkey
            statistics_dictkey = bitkey + '_' + conf_name_dict[i]
            statistics_dict[statistics_dictkey] = bitkey.upper(), conf_name_dict[i], noZ_stats(key)
            
        else:
            key = 'df_conf' + "".join(conf_name_dict[i].split("+")) + '_' + bitkey
            statistics_dictkey = bitkey + '_' + conf_name_dict[i]
            statistics_dict[statistics_dictkey] = bitkey.upper(), conf_name_dict[i], summary_stats(key)
			
#Print out the summary table
print('{:5s} {:8s} {:8s} {:15s} {:10s} {:15s} {:8s}'.format('Alg', \
                                                                    'Conf',\
                                                                'Nobj', '% in mask',\
                                                            '%in alg', '% high in alg', \
                                                            '% med in alg'))

print("-------------------------------------------------------------------------------")
print("-------------------------------------------------------------------------------")
for k, v in statistics_dict.items():
    alg, conf, stats = v
    nobj, permask, peralg, perhigh, permed = stats
    #print(type(permed))
    print('{:5s} {:8s} {:9s} {:15s} {:11s} {:15s} {:8s}'.format(alg, \
                                                                    conf,\
                                                                str(nobj), str(permask),\
                                                            str(peralg), str(perhigh), \
                                                                    str(permed)))
    if(k.split("_")[-1] == 'noZ'):
        print("-------------------------------------------------------------------------------")