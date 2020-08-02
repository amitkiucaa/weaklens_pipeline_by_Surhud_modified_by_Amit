import time
import sys
import os
import math
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d,interp2d
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
from scipy.integrate import dblquad
from astropy.io import fits
from scipy.integrate import quad, fixed_quad
import matplotlib.pyplot as plt
import matplotlib.pyplot as pltt
import emcee
from emcee.mpi_pool import MPIPool
import matplotlib.pyplot as pl
import corner
import sys
from time import sleep
import pandas as pd
import csv
import matplotlib.pyplot as pl
import corner
b_min=0.7# #bin min
b_max=0.9
dof=14


fffname="/mnt/home/student/camit/github/weaklens_pipeline/hsc_redmapper_overlap/output/Dsigma_%0.1f_%0.1f.dat"%(b_min,b_max)
sat_obs,sat_axis,sat_err,sat_pairs=np.loadtxt(fffname,usecols=(5,7,11,13),unpack=True)#z_oobs = observed delta sigma

f2ffname="/mnt/home/student/camit/github/weaklens_pipeline/redmapper/output/redmapper_satellites/run20/Dsigma_%0.1f_%0.1f.dat"%(b_min,b_max)
sat_obs2,sat_axis2,sat_err2=np.loadtxt(f2ffname,usecols=(5,7,11),unpack=True)
    

#writting covariance_matrix in a text file if dont exist
if os.path.exists("/mnt/home/student/camit/Data_Store/hsc_redmapper_overlap/covarinace_matrix_%0.1f_%0.1f.txt"%(b_min,b_max)):
    variance_matrix = np.loadtxt("/mnt/home/student/camit/Data_Store/hsc_redmapper_overlap/covarinace_matrix_%0.1f_%0.1f.txt"%(b_min,b_max), dtype='f', delimiter=' ')
else:
    #define covarince matrices
    #
    
    dfdict = {}
    for i in range(0, 320):                                     # we have 320 files of data
        ffname="/mnt/home/student/camit/github/weaklens_pipeline/hsc_redmapper_overlap/output/random_rotations/%0.1f_%0.1f/Dsigma.dat%05d"%(b_min,b_max,i)
        dfdict['%d' % i] = pd.read_csv(ffname, delim_whitespace=True, usecols=([5,7,11,14]), header=None, names=(["dsigma","radial_dist","error_data","R2_selection_bias"]),comment='#')
    
    num_lines = len(dfdict['%d' % 10].dsigma.values)      #number of lines in file for 10th  file.....10 is random number
    
    #mean delta sigma at for same radial bin in all files
    def mean_delta_sigma(rj):
        sum=0.0
        for ii in range(0,320):                 #division by 1.0+dfdict[] is done to apply R2 selection bias in the next line
            sum += dfdict['%d' % ii].dsigma.values[rj] / (1.0 + dfdict['%d' % ii].R2_selection_bias.values[rj])
    
        return sum/320.0
    
    #
    
    #variance
    def variance_j_k(j,k):                                     # j & k ranges 1 to num_lines
        mean_j=mean_delta_sigma(j)
        mean_k=mean_delta_sigma(k)
        sum_for_std_j_k=0.0
        for zz in range (0,320):                        #division by 1.0+dfdict[] is done to apply R2 selection bias in the next line
            sigma_diff_j = (dfdict['%d' % zz].dsigma.values[j] / (1.0 + dfdict['%d' % zz].R2_selection_bias.values[j])) - mean_j
            sigma_diff_k = (dfdict['%d' % zz].dsigma.values[k] / (1.0 + dfdict['%d' % zz].R2_selection_bias.values[k])) - mean_k
            sum_for_std_j_k += sigma_diff_j * sigma_diff_k
    
        return sum_for_std_j_k/320.0
    
    variance_matrix = []
    for iii in range(0,num_lines):  # A for loop for row entries
        a = []
        for jjj in range(0,num_lines):  # A for loop for column entries
            a.append(variance_j_k(iii,jjj))
        variance_matrix.append(a)
    np.savetxt('/mnt/home/student/camit/Data_Store/hsc_redmapper_overlap/covarinace_matrix_%0.1f_%0.1f.txt'%(b_min,b_max),np.matrix(variance_matrix),fmt='%.8f')

error_obs=np.zeros(len(variance_matrix[0]))
for dd in range(len(variance_matrix[0])):
    error_obs[dd]=np.sqrt(variance_matrix[dd][dd])

#print error_obs

variance_matrix2 = np.loadtxt("/mnt/home/student/camit/Data_Store/covarinace_matrix_%0.1f_%0.1f.txt"%(b_min,b_max), dtype='f', delimiter=' ')
error_obs2=np.zeros(len(variance_matrix2[0]))
for dd in range(len(variance_matrix2[0])):
    error_obs2[dd]=np.sqrt(variance_matrix2[dd][dd])


print (error_obs2,"\n",error_obs,"\n",sat_obs2,"\n",sat_obs)




