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
np.random.seed(10)
os.chdir(r"/mnt/home/student/camit/github/HSCWLsimulations/dark_emulator") 

#

b_min=0.7# #bin min
b_max=0.9
dof=14
Nwalkers,Niter_burn,Niter=300,5000,15000  
#Nwalkers,Niter_burn,Niter=50,20,20  
#m_parent_in,m_sat_in,m_stellar_in,conc_parent_in,conc_sat_in,mck_in=14.303,  12.059,  10.685,   4.466,   1.771,   0.019          #0.1-0.3
#m_parent_in,m_sat_in,m_stellar_in,conc_parent_in,conc_sat_in,mck_in=14.317,  12.365,  10.258,   3.819,   5.205,   0.033          #0.3-0.5
#m_parent_in,m_sat_in,m_stellar_in,conc_parent_in,conc_sat_in,mck_in=14.358,  11.963,  7.711,    3.484,   9.134,   0.035          #0.5-0.7
m_parent_in,m_sat_in,m_stellar_in,conc_parent_in,conc_sat_in,mck_in=14.37,   12.233,  6.281,    3.903,   13.052,  0.127          #0.7-0.9

#

#parameters from mcmc another mcmc sampling done for field galaxies
field_galaxies_mass=10**12.91
field_galaxies_stellar_mass=10**9.94
field_galaxies_conc=3.39
#Hand plot parameters
config_rmin=-2  #in logspace  from weaklens pipeline
config_rmax=np.log10(5)#in logspace
config_bin=20

#
#m_parent_plot,m_sat_plot,m_stellar_plot,conc_parent_plot,conc_sat_plot,mck_plot=14.303,  12.059,  10.685,   4.466,   1.771,   0.019    #0.1-0.3
#m_parent_plot,m_sat_plot,m_stellar_plot,conc_parent_plot,conc_sat_plot,mck_plot= 14.317, 12.365, 10.258,  3.819,  5.205,  0.033       #0.3-0.5
#m_parent_plot,m_sat_plot,m_stellar_plot,conc_parent_plot,conc_sat_plot,mck_plot= 14.358, 11.963,  7.711,  3.484,  9.134,  0.035#0.5-0.7
m_parent_plot,m_sat_plot,m_stellar_plot,conc_parent_plot,conc_sat_plot,mck_plot= 14.37,  12.233,  6.281,  3.903,  13.052,  0.127#0.7-0.9


theta_hand=m_parent_plot,m_sat_plot,m_stellar_plot,conc_parent_plot,conc_sat_plot,mck_plot

fffname="/mnt/home/student/camit/github/weaklens_pipeline/hsc_redmapper_overlap/output/Dsigma_%0.1f_%0.1f.dat"%(b_min,b_max)
sat_obs,sat_axis,sat_err,sat_pairs=np.loadtxt(fffname,usecols=(5,7,11,13),unpack=True)#z_oobs = observed delta sigma

f2ffname="/mnt/home/student/camit/github/weaklens_pipeline/redmapper/output/redmapper_satellites/run20/Dsigma_%0.1f_%0.1f.dat"%(b_min,b_max)
sat_obs2,sat_axis2,sat_err2=np.loadtxt(f2ffname,usecols=(5,7,11),unpack=True)
    

f1ffname="/mnt/home/student/camit/github/weaklens_pipeline/redmapper/output/redmapper_satellites/run20/stellar_counterparts/Dsigma_%0.1f_%0.1f.dat"%(b_min,b_max)
stellar_obs,stellar_axis,stellar_err=np.loadtxt(f1ffname,usecols=(5,7,11),unpack=True)
    
pairs_plot=np.logspace(config_rmin,config_rmax,config_bin+1) #plotting for weighted average taking pairs into consideration #+1 so that r1+r2/2 will give same ans
spline_totalpairs=interp1d(np.log10(sat_axis),np.log10(sat_pairs/np.diff(pairs_plot)), fill_value="extrapolate")

#

x=np.loadtxt("/mnt/home/student/camit/Data_Store/required_text_for_storing_galaxies_distance_0.05_1.2_from_center.txt",usecols=0,unpack=True)
n, bin, patches = plt.hist(x, bins =100,histtype="step",density=True)
bin_center = bin[:-1] + np.diff(bin) / 2
number_density_spline=interp1d(bin_center,n,kind="cubic")
factor=fixed_quad(lambda Rsat: number_density_spline(Rsat), b_min, b_max,n=100)[0]
plt.close()  #otherwise plot of number desity will be also there along wity delta sigma plot which we are going ot plot


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

cov=variance_matrix
icov = np.linalg.inv(cov)
#cov = np.triu(cov)
#cov += cov.T - np.diag(cov.diagonal())
#cov = np.dot(cov, cov)

if os.path.exists("/mnt/home/student/camit/Data_Store/tcovarinace_matrix_%0.1f_%0.1f.txt"%(b_min,b_max)):
    tvariance_matrix = np.loadtxt("/mnt/home/student/camit/Data_Store/tcovarinace_matrix_%0.1f_%0.1f.txt"%(b_min,b_max), dtype='f', delimiter=' ')
else:

    tdfdict = {}                                                 
    for i in range(0, 320):                                     # we have 320 files of data
        tffname="/mnt/home/student/camit/github/weaklens_pipeline/redmapper/output/redmapper_satellites/run20/stellar_counterparts/random_rotations/%0.1f_%0.1f/Dsigma.dat%05d"%(b_min,b_max,i)
        tdfdict['%d' % i] = pd.read_csv(tffname, delim_whitespace=True, usecols=([5,7,11,14]), header=None, names=(["tdsigma","tradial_dist","terror_data","tR2_selection_bias"]),comment='#')
    
    #
    
    tnum_lines = len(tdfdict['%d' % 10].tdsigma.values)      #number of lines in file for 10th  file.....10 is random number
    
    #
    
    
    #mean delta sigma at for same radial bin in all files
    def tmean_delta_sigma(rj):
        tsum=0.0
        for ii in range(0,320):                 #division by 1.0+dfdict[] is done to apply R2 selection bias in the next line
            tsum += tdfdict['%d' % ii].tdsigma.values[rj] / (1.0 + tdfdict['%d' % ii].tR2_selection_bias.values[rj])
    
        return tsum/320.0
    
    #
    
    #variance
    def tvariance_j_k(j,k):                                     # j & k ranges 1 to num_lines
        tmean_j=tmean_delta_sigma(j)
        tmean_k=tmean_delta_sigma(k)
        tsum_for_std_j_k=0.0
        for zz in range (0,320):                        #division by 1.0+dfdict[] is done to apply R2 selection bias in the next line
            tsigma_diff_j = (tdfdict['%d' % zz].tdsigma.values[j] / (1.0 + tdfdict['%d' % zz].tR2_selection_bias.values[j])) - tmean_j
            tsigma_diff_k = (tdfdict['%d' % zz].tdsigma.values[k] / (1.0 + tdfdict['%d' % zz].tR2_selection_bias.values[k])) - tmean_k
            tsum_for_std_j_k += tsigma_diff_j * tsigma_diff_k
        return tsum_for_std_j_k/320.0
    
    #
    
    #code for matrix of variance
    tvariance_matrix = []
    for iii in range(0,tnum_lines):  # A for loop for row entries
        ta = []
        for jjj in range(0,tnum_lines):  # A for loop for column entries
            ta.append(tvariance_j_k(iii,jjj))
        tvariance_matrix.append(ta)
    np.savetxt('/mnt/home/student/camit/Data_Store/tcovarinace_matrix_%0.1f_%0.1f.txt'%(b_min,b_max),np.matrix(tvariance_matrix),fmt='%.8f')

terror_obs=np.zeros(len(tvariance_matrix[0]))
for dd in range(len(tvariance_matrix[0])):
    terror_obs[dd]=np.sqrt(tvariance_matrix[dd][dd])
#tcov=tvariance_matrix
#ticov = np.linalg.inv(tcov)

#

class DeltaSigmaSatellite:    #    nfw analytical profile #give credit to the paper
    def __init__(self,halo_mass,conc):   # mass and concentration are the parameters that difine a nfw profie
        self.M200=halo_mass
        self.c200=conc 
  
        self.Omm = 0.315
        self.gee = 4.3e-9
        self.rhocrit = 3e4/(8.0*np.pi*self.gee)
        self.rhobar = self.rhocrit*self.Omm

        self.R200 = (self.M200/(4./3.*np.pi*self.rhobar*200))**(1./3.)
        self.rs = self.R200/self.c200

        def mu(x):
            return np.log(1.+x)-x/(1.+x)

        self.rhos = self.M200/(4.0*np.pi*self.rs**3*mu(self.c200))


    def g(self,x):
        res = x*1.0

        idx = x<1
        xlt1 = x[idx]
        fac = (1-xlt1)/(1+xlt1)
        res[idx] = 8*np.arctanh(fac**0.5)/(xlt1**2*(1-xlt1**2)**0.5) + 4/xlt1**2*np.log(xlt1/2) - 2/(xlt1**2-1) + 4*np.arctanh(fac**0.5)/(xlt1**2-1)/(1-xlt1**2)**0.5

        idx = x>1
        xgt1 = x[idx]
        fac = (xgt1-1)/(1+xgt1)
        res[idx] = 8*np.arctan(fac**0.5)/(xgt1**2*(xgt1**2-1)**0.5) + 4/xgt1**2*np.log(xgt1/2) - 2/(xgt1**2-1) + 4*np.arctan(fac**0.5)/(xgt1**2-1)**1.5

        idx= x==1
        res[idx]=10./3. + 4*np.log(0.5)

        return res
    
    
    def dsigma(self,r):
        x = r/self.rs

        return self.rs*self.rhos*self.g(x)/1e12

    

#for Host    
   

class DeltaSigmaParent_initial:
    def __init__(self,Mparent,parent_conc,mck):  # mass,conc for nfw & mck for miscentering_kernal
        self.Mparent=Mparent
        self.conc=parent_conc
        self.nfw_parent=DeltaSigmaSatellite(self.Mparent,self.conc)
        self.rmin_spline=-1.5
        self.initialized=False
        self.mck=mck
        rs = np.logspace(self.rmin_spline,1.1,50) #sampling Range
        self.delta_sigma_parent_spline=interp1d(rs, self.nfw_parent.dsigma(rs), kind="cubic")


    def Delta_sigma_parent(self,Rp,RSAT):  #take Rp as a number, RSAT can be an array or number   # By our method
        def integrand(thetat,Rsat):  #integrate delta sigma(r) on theta 0 to np.pi
            if np.isscalar(thetat):
                theta = np.array([thetat])
            else:

                theta = thetat

            cost=np.cos(theta)
            r=np.sqrt(Rp**2 + Rsat**2 + 2*Rp*Rsat*cost )
            cosalpha=(Rp+(Rsat*cost)) / r
            ans = np.zeros(theta.size)   
            idx = np.log10(r)<self.rmin_spline
            ans[idx] = self.delta_sigma_parent_spline(10**self.rmin_spline)*(2*(cosalpha[idx]**2)-1)
            ans[~idx] = self.delta_sigma_parent_spline(r[~idx])*(2*(cosalpha[~idx]**2)-1)

            if np.isscalar(thetat):
                return ans[0]
            else:
                return ans

        if np.isscalar(RSAT):
            Rsatt = np.array([RSAT])
        else:
            Rsatt = RSAT

        integration_fquad=np.zeros(Rsatt.size)

        for i in range(Rsatt.size):
            integration_fquad[i]= fixed_quad(integrand,0,np.pi, n=50,args=[Rsatt[i]])[0]

        if np.isscalar(RSAT):
            return integration_fquad[0]/(np.pi)

        return integration_fquad/(np.pi)

class DeltaSigmaParent():
    def __init__(self,Mparent,parent_conc,mck):  # mass,conc for nfw & mck for miscentering_kernal
        self.Mparent=Mparent
        self.conc=parent_conc
        self.nfw_parent=DeltaSigmaSatellite(self.Mparent,self.conc)
        self.rmin_spline=-1.5
        self.initialized=False
        self.mck=mck
        self.parent_class=DeltaSigmaParent_initial(self.Mparent,self.conc,self.mck)

    def miscentered_probablity(self,Rsat,Rsat0,sigma):
        return (1/(2*np.pi*Rsat0))*(np.exp(-(Rsat-Rsat0)**2/(2*sigma**2))-np.exp(-(Rsat+Rsat0)**2/(2*sigma**2)) )

    def dsigma_parent_miscentered_surhud(self,R,Rsat,sigma):
    
        #sampling_axis=np.linspace(slim_of_integration_min,slim_of_integration_max,50)
        sampling_axis=np.linspace(0.001,2*b_max,100)
        Delta_sigma_parent_spline=interp1d(sampling_axis,self.parent_class.Delta_sigma_parent(R,sampling_axis))
        if np.isscalar(Rsat):
            RRsat=np.array([Rsat])
        else:
            RRsat=Rsat
        prob_norm_const=np.zeros(RRsat.size)
        ans=np.zeros(RRsat.size)
 

        def integrand(r,r0,sigma1):
            return Delta_sigma_parent_spline(r)*self.miscentered_probablity(r,r0,sigma1)
       
        for i in range(RRsat.size):
             p=3*sigma
             rsat=RRsat[i]
             lim_of_integration_min=rsat-p
             lim_of_integration_max=rsat+p
             ans[i]=fixed_quad(integrand,lim_of_integration_min,lim_of_integration_max,args=[rsat,sigma],n=30)[0]
             prob_norm_const[i]= fixed_quad(self.miscentered_probablity,0.01,lim_of_integration_max,args=[rsat,sigma],n=30)[0]
        if np.isscalar(Rsat):
            return ans[0]/prob_norm_const[0]
        return ans/prob_norm_const



    def weighted_Delta_sigma_fixed_parent_final_miscentered_surhud(self,RR):            
        def integrand(Rsat,R,sigma):
            return number_density_spline(Rsat)*self.dsigma_parent_miscentered_surhud(R,Rsat,sigma)
        if np.isscalar(RR):
            RRR = np.array([RR])
        else:
            RRR = RR       

        ans =np.zeros(RRR.size)
        for i in range(RRR.size):
            ans[i]=fixed_quad(integrand,b_min,b_max,args=[RRR[i],self.mck],n=50)[0]   
        if np.isscalar(RR):
            return ans[0]/factor
        return ans/factor 



    def weighted_Delta_sigma_fixed_parent_final(self,RR):
        def integrand(Rsat,R):
            return number_density_spline(Rsat)*self.Delta_sigma_parent(R,Rsat)
        if np.isscalar(RR):
            RRR = np.array([RR])
        else:
            RRR = RR

        ans =np.zeros(RRR.size)
        for i in range(RRR.size):
            ans[i]=fixed_quad(integrand,b_min,b_max,args=[RRR[i]])[0]
        return ans/factor




    def dsigma_miscentered(self,rrr):
        return self.weighted_Delta_sigma_fixed_parent_final_miscentered_surhud(rrr)    

    def dsigma(self,rrr):
        return self.weighted_Delta_sigma_fixed_parent_final(rrr)
  

#for stellar    
    
class DeltaSigmaStellar:
    def __init__(self,Mstellar):
        self.Mstellar=Mstellar
    
    def dsigma(self,R):
        if np.isscalar(R):
            rrrr=np.array([R])
        else:
            rrrr=R
            
        return (self.Mstellar/(np.pi*(rrrr**2)))/1e12    



    
class mcmc_pool:
    def __init__(self,bin_min,bin_max):
        self.obs=sat_obs#observatinal data
        self.Nobs = self.obs.size #size of array
        self.obs_axis=sat_axis#np.logspace(-2,np.log10(5),self.Nobs) #plotting range
        self.rbin_min=bin_min
        self.rbin_max=bin_max
        self.No_of_observed_points=len(self.obs)    

    #theoretical model

    def theoretical_prediction_mcmc(self,RRp,theta,rbin_min,rbin_max):    #model for mcmc
        lnmparent_true,lnmsat_true,lnmstellar_true,conc_parent,conc_sat,mck=theta
        dsigma_sat=DeltaSigmaSatellite(10**lnmsat_true,conc_sat).dsigma(RRp)
        dsigma_parent=DeltaSigmaParent(10**lnmparent_true,conc_parent,mck).dsigma_miscentered(RRp)#mck ->miscentering kernal
        dsigma_stellar=DeltaSigmaStellar(10**lnmstellar_true).dsigma(RRp)
        return dsigma_sat+dsigma_stellar+dsigma_parent


    def theoretical_prediction_mcmc9(self,RRp,theta,rbin_min,rbin_max):    #model for mcmc
        lnmparent_true,lnmsat_true,lnmstellar_true,conc_parent,conc_sat,mck=theta
        dsigma_sat=DeltaSigmaSatellite(10**lnmsat_true,conc_sat).dsigma(RRp)
        dsigma_parent=DeltaSigmaParent(10**lnmparent_true,conc_parent,mck).dsigma_miscentered(RRp)#mck ->miscentering kernal
        dsigma_stellar=DeltaSigmaStellar(10**lnmstellar_true).dsigma(RRp)
        return dsigma_sat,dsigma_stellar,dsigma_parent




     

    def theoretical_prediction1(self,RRp,theta,rbin_min,rbin_max):    #miscentering model
        lnmparent_true,lnmsat_true,lnmstellar_true,conc_parent,conc_sat,mck=theta
        dsigma_sat=DeltaSigmaSatellite(10**lnmsat_true,conc_sat).dsigma(RRp)
        dsigma_parent=DeltaSigmaParent(10**lnmparent_true,conc_parent,mck).dsigma_miscentered(RRp)#mck ->miscentering kernal
        dsigma_stellar=DeltaSigmaStellar(10**lnmstellar_true).dsigma(RRp)
        return dsigma_parent+dsigma_sat+dsigma_stellar
        #return dsigma_parent



    def theoretical_prediction2(self,RRp,theta,rbin_min,rbin_max):    #without miscentering 
        lnmparent_true,lnmsat_true,lnmstellar_true,conc_parent,conc_sat,mck=theta
        dsigma_sat=DeltaSigmaSatellite(10**lnmsat_true,conc_sat).dsigma(RRp)
        dsigma_parent=DeltaSigmaParent(10**lnmparent_true,conc_parent,mck).dsigma(RRp)  #no miscentering
        dsigma_stellar=DeltaSigmaStellar(10**lnmstellar_true).dsigma(RRp)
        return dsigma_parent+dsigma_sat+dsigma_stellar
        #return dsigma_parent


    def theoretical_prediction3(self,RRp,theta,rbin_min,rbin_max):    #miscentering + cen-sat_pairs model
        lnmparent_true,lnmsat_true,lnmstellar_true,conc_parent,conc_sat,mck=theta
        dsigma_sat=DeltaSigmaSatellite(10**lnmsat_true,conc_sat).dsigma(RRp)
        dsigma_parent=DeltaSigmaParent(10**lnmparent_true,conc_parent,mck).dsigma_miscentered(RRp)#mck ->miscentering kernal
        dsigma_stellar=DeltaSigmaStellar(10**lnmstellar_true).dsigma(RRp)
        
        dsigma_thp=dsigma_parent+dsigma_sat+dsigma_stellar
        spline_thpredict=interp1d(np.log10(RRp),dsigma_thp, fill_value="extrapolate")
        
        storing_delta_sigma1=np.zeros(len(pairs_plot)-1)
        if not len(RRp)==len(storing_delta_sigma1):
            print("You are not plotting with observational axis")
            sys.exit()
        
        for i in range(len(pairs_plot)-1):
            a=pairs_plot[i]
            b=pairs_plot[i+1]
            integration_numerator=fixed_quad(lambda x: (10**spline_totalpairs(np.log10(x)))*spline_thpredict(np.log10(x)),a,b,n=50)[0]
            integration_denominator=fixed_quad(lambda x: (10**spline_totalpairs(np.log10(x))),a,b,n=50)[0]
            delta_sigma_a_b=integration_numerator/integration_denominator
            storing_delta_sigma1[i] =delta_sigma_a_b
            #integration_numerator2=quad(lambda x: spline_thpredict(np.log10(x)),a,b)[0]
            #integration_denominator2=b-a
            #delta_sigma_a_b2=integration_numerator2/integration_denominator2
            #storing_delta_sigma2[i] =delta_sigma_a_b2
            #storing_sat_axis[i]=(a+b)/2
            #print("delta_sigma averaged,not averaged,observed",delta_sigma_a_b,delta_sigma_a_b2)
        return storing_delta_sigma1


    def theoretical_prediction4(self,RRp,theta,rbin_min,rbin_max):    # cen-sat_pairs model, no miscentering
        lnmparent_true,lnmsat_true,lnmstellar_true,conc_parent,conc_sat,mck=theta
        dsigma_sat=DeltaSigmaSatellite(10**lnmsat_true,conc_sat).dsigma(RRp)
        dsigma_parent=DeltaSigmaParent(10**lnmparent_true,conc_parent).dsigma(RRp)
        dsigma_stellar=DeltaSigmaStellar(10**lnmstellar_true).dsigma(RRp)
        
        dsigma_thp=dsigma_parent+dsigma_sat+dsigma_stellar
        spline_thpredict=interp1d(np.log10(RRp),dsigma_thp, fill_value="extrapolate")

        
        storing_delta_sigma1=np.zeros(len(pairs_plot)-1)
        if not len(RRp)==len(storing_delta_sigma1):
            print("You are not plotting with observational axis")
            sys.exit()
        

        for i in range(len(pairs_plot)-1):
            a=pairs_plot[i]
            b=pairs_plot[i+1]
            integration_numerator=fixed_quad(lambda x: (10**spline_totalpairs(np.log10(x)))*spline_thpredict(np.log10(x)),a,b,n=50)[0]
            integration_denominator=fixed_quad(lambda x: (10**spline_totalpairs(np.log10(x))),a,b,n=50)[0]
            delta_sigma_a_b=integration_numerator/integration_denominator
            storing_delta_sigma1[i] =delta_sigma_a_b
            #integration_numerator2=quad(lambda x: spline_thpredict(np.log10(x)),a,b)[0]
            #integration_denominator2=b-a
            #delta_sigma_a_b2=integration_numerator2/integration_denominator2
            #storing_delta_sigma2[i] =delta_sigma_a_b2
            #storing_sat_axis[i]=(a+b)/2
            #print("delta_sigma averaged,not averaged,observed",delta_sigma_a_b,delta_sigma_a_b2)
        return storing_delta_sigma1
    



    def lnprior(self,theta):
        lnmparent_t,lnmsat_t,lnmstellar_t,conc_parent_t,conc_sat_t,mck_t = theta    #parameters gets values in p sequntly
    
        if  (10<= lnmparent_t < 16) and (8<= lnmsat_t < 16) and (5<=lnmstellar_t < 14)and (1<=conc_parent_t < 40)and (1<=conc_sat_t < 40)and (0.001<mck_t<b_min/3.1):
            return 0.0
        return -np.inf

    #posterior
    def lnprob(self,theta):
        lnmparent_t,lnmsat_t,lnmstellar_t,conc_parent_t,conc_sat_t,mck_t = theta
        

        lp = self.lnprior(theta)
         
        if not np.isfinite(lp):
            return -np.inf,theta,5.0,np.ones(self.No_of_observed_points) #0 for r_chisquare
        
        no_of_free_parameters=len(theta)
        dof=self.No_of_observed_points-no_of_free_parameters 
        
        model=self.theoretical_prediction_mcmc(self.obs_axis,theta,self.rbin_min,self.rbin_max)
        
        diff=(self.obs - model)
        aaaa= np.dot(icov,diff)   
        chisq = np.dot(diff,aaaa)   #np.dot(Delta, np.dot(icov, Delta)) 
        ln_prior= -0.5*chisq  
        if not np.isfinite(chisq):
            print("encounterd nan or infintie for lnprior",model,theta,aaaa)
            return -np.inf,theta,5.0,np.ones(self.No_of_observed_points) #0 for chisquare
        #print ("ln_Mparent,ln_Msatellite,ln_Mstellar,conc_parent, conc_sat, miscentering_kernal Chisq/DOF  :  %f \t%f \t %f \t%f \t%f \t%f  \t%f "%(theta[0],theta[1],theta[2],theta[3],theta[4],theta[5],chisq/dof))
        print ("%f \t%f \t %f \t%f \t%f \t%f  \t%f "%(theta[0],theta[1],theta[2],theta[3],theta[4],theta[5],chisq/dof))
        #print(theta,chisq/dof)
        return lp+ln_prior,theta,chisq/dof,model


bb=mcmc_pool(b_min,b_max)
def surhudprob(theta):
    return bb.lnprob(theta)  


if __name__ == "__main__":
        
    obs_axis=sat_axis
    obs=sat_obs

    amit=mcmc_pool(b_min,b_max)
    #plotting for final figure theta_max by hand
    mod=amit.theoretical_prediction_mcmc(sat_axis,theta_hand,b_min,b_max)
    mod_field=DeltaSigmaSatellite(field_galaxies_mass,field_galaxies_conc).dsigma(sat_axis)+ DeltaSigmaStellar(field_galaxies_stellar_mass).dsigma(sat_axis)



    plt.plot(sat_axis,mod,"-",label="Satellite_galaxies_Model")
    plt.errorbar(sat_axis,mod_field,fmt=":",label="Field_galaxies_Model")
    plt.errorbar(sat_axis,sat_obs,yerr=error_obs,fmt="*",mfc="grey",label="Satellite_Galaxies_Obs")     
    plt.errorbar(stellar_axis,stellar_obs,yerr=terror_obs,fmt="p",mfc="grey",ms=4,label="Field Gallaxies_Obs")
    plt.errorbar(sat_axis2,sat_obs2,yerr=sat_err2,fmt="p",mfc="black",ms=4,label="sat old")
    #plt.plot(sat_axis,mod-mod_field,label="enviorenmental_effects")
    plt.legend()
    plt.xscale('log')
    #plt.yscale('log')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
        
    plt.ylabel(r"$\overline{\Delta\Sigma}$ [hM$_\odot$/pc$^2$]",fontsize=15)
    plt.xlabel(r"$R [h^{-1} Mpc$]",fontsize=15)
    plt.legend(loc="upper right",prop={'size':13})
    #plt.ylim(-50, 180) 
    #plt.axhline(0,ls="--")
    plt.title("%0.1f_%0.1f"%(b_min,b_max),size='x-large')
        
 




           
    diff=(mod - sat_obs)
    aaaa= np.dot(icov,diff)   
    chisq_corr = np.dot(diff,aaaa)#
    plt.text(1,-35,"$\chi^2_{red}$ = %0.2f"%(chisq_corr/dof), bbox=dict(facecolor='grey', alpha=0.2),fontsize=15)




    """axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
    axins.imshow(Z2, extent=extent, interpolation="nearest",origin="lower")
    # sub region of the original image
    x1, x2, y1, y2 = -1.5, -0.9, -2.5, -1.9
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    axins.set_xticklabels('')
    axins.set_yticklabels('')

    ax.indicate_inset_zoom(axins)"""





    #plt.plot(sat_axis,(mod-sat_obs)/error_obs,"*",label="chi")
    #plt.savefig("/mnt/home/student/camit/outputs/maximum_probable_model%0.1f_%0.1f_testingg.png"%(b_min,b_max))

    #plt.savefig("/mnt/home/student/camit/outputs/%0.1f_%0.1f_testing.png"%(b_min,b_max))
    plt.savefig("/mnt/home/student/camit/outputs/maximum_probable_model%0.1f_%0.1f_testingg.png"%(b_min,b_max))
    #plt.savefig("/mnt/home/student/camit/outputs_sat_conc_fixed/maximum_probable_model%0.1f_%0.1f_testingg.png"%(b_min,b_max))
    plt.savefig("maximum_probable_model%0.1f_%0.1f_testingg.png"%(b_min,b_max))   
    plt.show() 
    #plt.clf()
    #plt.close()
    print("chisq is", chisq_corr/dof)


    """

    plt.close()


    mod=amit.theoretical_prediction_mcmc(sat_axis,theta_hand,b_min,b_max)


    
    plt.plot(sat_axis,mod,"-",label="Satellite_galaxies_Model")
    plt.errorbar(sat_axis,mod_field,fmt=":",label="Field_galaxies_Model")
    plt.errorbar(sat_axis,sat_obs,yerr=error_obs,fmt="*",label="Satellite_Galaxies_Obs")     
    plt.errorbar(stellar_axis,stellar_obs,yerr=terror_obs,fmt="p",label="Field Gallaxies_Obs")
    #plt.plot(sat_axis,mod-mod_field,label="enviorenmental_effects")
    plt.legend()
    plt.xscale('log')
    #plt.yscale('log')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
        
    plt.ylabel(r"$\overline{\Delta\Sigma}$ [hM$_\odot$/pc$^2$]",fontsize=15)
    plt.xlabel(r"$R [h^{-1} Mpc$]",fontsize=15)
    plt.legend(loc="upper right",prop={'size':13})
    plt.ylim(-50, 180) 
    #plt.axhline(0,ls="--")
    plt.title("Model_vs_Observations_for_galaxies_%0.1f_%0.1f_$h^{-1}$ Mpc_from_cluster_center"%(b_min,b_max),size='x-large')
    """    
            




    print("Done")
    

 
