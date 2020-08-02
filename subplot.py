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
np.random.seed(10)

#
fffname="/mnt/home/student/camit/github/weaklens_pipeline/arindam/output/Dsigma.dat"
sat_obs,sat_axis,sat_err,sat_pairs=np.loadtxt(fffname,usecols=(5,7,11,13),unpack=True)#z_oobs = observed delta sigma

#writting covariance_matrix in a text file if dont exist
if os.path.exists("/mnt/home/student/camit/Data_Store/covarinace_matrix_arindam_lenses.txt"):
    variance_matrix = np.loadtxt("/mnt/home/student/camit/Data_Store/covarinace_matrix_arindam_lenses.txt", dtype='f', delimiter=' ')
else:
    #define covarince matrices
    #
    print("Processing")    
    dfdict = {}
    for i in range(0, 320):                                     # we have 320 files of data
        ffname="/mnt/home/student/camit/github/weaklens_pipeline/arindam/output/covarinace_matrices/Dsigma.dat%05d"%(i)
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
    np.savetxt("/mnt/home/student/camit/Data_Store/covarinace_matrix_arindam_lenses.txt",np.matrix(variance_matrix),fmt='%.8f')

error_obs=np.zeros(len(variance_matrix[0]))
for dd in range(len(variance_matrix[0])):
    error_obs[dd]=np.sqrt(variance_matrix[dd][dd])


plt.figure(figsize=(10,5))

ax=plt.subplot(1,2,1)
ax.errorbar(sat_axis,sat_obs,yerr=error_obs,fmt="*",ls="-",label="loglog with errorbars")     
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_ylabel(r"$\overline{\Delta\Sigma}$ [hM$_\odot$/pc$^2$]",fontsize=15)
ax.set_xlabel(r"$R [h^{-1} Mpc$]",fontsize=15)
ax.legend(loc="upper right",prop={'size':13})



ax=plt.subplot(1,2,2)
ax.errorbar(sat_axis,sat_obs,yerr=error_obs,fmt="*",ls="-",label="semilog with errorbars")     
ax.set_xscale('log')

#ax.set_ylabel(r"$\overline{\Delta\Sigma}$ [hM$_\odot$/pc$^2$]",fontsize=15)
ax.set_xlabel(r"$R [h^{-1} Mpc$]",fontsize=15)
ax.legend(loc="upper right",prop={'size':13})


"""
ax=plt.subplot(2,2,3)
ax.errorbar(sat_axis,sat_obs,fmt="*",ls="-",label="loglog without errorbars")     
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_ylabel(r"$\overline{\Delta\Sigma}$ [hM$_\odot$/pc$^2$]",fontsize=15)
ax.set_xlabel(r"$R [h^{-1} Mpc$]",fontsize=15)
ax.legend(loc="upper right",prop={'size':13})



ax=plt.subplot(2,2,4)
#ax.errorbar(sat_axis,sat_obs,fmt="*",label="semilog without errorbars")     
ax.set_xscale('log')

#ax.set_ylabel(r"$\overline{\Delta\Sigma}$ [hM$_\odot$/pc$^2$]",fontsize=15)
#ax.legend(loc="upper right",prop={'size':13})
ax.set_xlabel(r"$R [h^{-1} Mpc$]",fontsize=15)
"""

plt.savefig("subplot.png")
plt.show()

#plt.savefig("/mnt/home/student/camit/outputs/%0.1f_%0.1f_testing.png"%(b_min,b_max))






