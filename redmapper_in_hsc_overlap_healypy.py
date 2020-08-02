import numpy as np	
from astropy.io import fits as fits
import healpy as hp
import time,os,sys
from tqdm import tqdm

#function for converting RA,DEC to Healypy Map index
def DecRaToIndex(NsidE,declination,rightassension):    
    return hp.pixelfunc.ang2pix(NsidE,np.radians(-declination+90.),np.radians(rightassension))

def IndexToDeclRa(Nside,index):
    theta,phi=hp.pixelfunc.pix2ang(Nside,index)
    return -np.degrees(theta-np.pi/2.),np.degrees(phi)


#reading HSC S16 map
print("Reading HSC Map")
Hsc_s16_map = hp.read_map("DataStore/S16A_fdfc_hp.fits").astype(np.bool)
# get Nside from this Map 
print("Done Reading HSC Map")


#to view map 
"""hp.mollview(
    Hsc_s16_map,
    coord=["G", "E"],
    title="Histogram equalized Ecliptic",
    unit="mK",
    norm="hist",
    min=-1,
    max=1,
)

hp.graticule()   #will add meridian lines"""


#import Redmapper catalog  RA and DEC


print("Importing Redmapper Data")
hdulist=fits.open("DataStore/redmapper/redmapper_member_modified.fits")
hdulist2 =fits.open("DataStore/redmapper/redmapper_dr8_public_v6.3_catalog.fits")
data = hdulist[1].data
data2 = hdulist2[1].data

#choose clusers where culster centre is defined with probablity more than 95%
tp_cen=data2["P_CEN"].astype("float64")
aa=np.array([i[0] for i in tp_cen])
tidx2=aa>0.95
data2=data2[tidx2]
cid_idx=data2["id"]  #in cluster catalog
mid_idx=data["id"]   #member catalog 
tidx = np.where(np.in1d(mid_idx, cid_idx),True,False)
data_members=data[tidx]



Ra=data_members["RA"]    # RA and DEC are in Degree
Dec=data_members["DEC"]
z=data_members["red_zcluster"]
#z=data_members["z"]
dist=data_members["R"]
#Pick Nside from above Map
Nside=hp.pixelfunc.get_nside(Hsc_s16_map)

print("Redmapper data loaded")


#Now creating Map for these RA Dec of Redmapper 
pixel_indices =DecRaToIndex(Nside,Dec,Ra)    #pixel indices for Redmapper
map_redmapper=np.zeros(hp.nside2npix(Nside))     #No of pixels corresponding to Nside
map_redmapper[np.unique(pixel_indices)]=1                   #writing True for pixels where Object is present



red_hsc_overlap = np.where((Hsc_s16_map==map_redmapper) & (Hsc_s16_map==1),True,False)


#to view map 
"""hp.mollview(
    red_hsc_overlap,
    coord=["G", "E"],
    title="HSC_Redmapper_overlap",
    unit="mK",
    norm="hist",
    min=-1,
    max=1,
)

hp.graticule()   #will add meridian lines"""

#creating index for red_hsc_overlap

index_hsc_over = np.where(red_hsc_overlap == True)
indx_red_hsc_overlap=index_hsc_over[0]

#writing Ra,Dec,Z if it lie in red_hsc_overlap
#making dictionary of RA only if they fall in hsc field


print("Writing overlap info to hsc_red_overlap.txt")

fp=open("hsc_red_members_overlap.txt","w")
mydict_new={}
for i in tqdm(range(Ra.size)):
    x=DecRaToIndex(Nside,Dec[i],Ra[i])
    if x in indx_red_hsc_overlap:
        mydict_new["%0.5f"%Ra[i]+"+"+"%0.5f"%Dec[i]+"+"+"%0.4f"%z[i]]=x
        fp.write("%0.5f,%0.5f,%0.4f,%0.4f\n"%(Ra[i],Dec[i],z[i],dist[i]))

fp.close()

"""
fp=open('state.txt','w')
for key in mydict_new.keys():
    a,b,c=key.split("+")
    fp.write("%0.5f,%0.5f,%0.4f\n"%(float(a),float(b),float(c)))
"""

print("\n**************\n    Done\n****************")


