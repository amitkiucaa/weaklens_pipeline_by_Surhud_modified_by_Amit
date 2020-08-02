from astropy.io import fits
import numpy as np

a,b,c,d=np.loadtxt("hsc_red_members_overlap.txt",usecols=(0,1,2,3),delimiter=",",unpack=True)

a1 = a#np.array(['NGC1001', 'NGC1002', 'NGC1003'])   #first column
a2 = b#np.array([11.1, 12.3, 15.2])                  #second column
a3=c
a4=d

col1 = fits.Column(name='redmapper_Ra', format='E', array=a1)
col2 = fits.Column(name='redmapper_Dec', format='E', array=a2)
col3 = fits.Column(name='redmapper_redshift', format='E', array=a3)
col4 = fits.Column(name='dist_from_cluster_center', format='E', array=a4)

cols = fits.ColDefs([col1, col2,col3,col4])
tbhdu = fits.BinTableHDU.from_columns(cols)

prihdr = fits.Header()

prihdr['Table_creater'] = 'Amit'
prihdr['COMMENT'] = "create hsc_redmapper_overlap fits file"

prihdu = fits.PrimaryHDU(header=prihdr)

thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto('redmapper_members_hsc_overlap_.fits')

print ("fits file created")
