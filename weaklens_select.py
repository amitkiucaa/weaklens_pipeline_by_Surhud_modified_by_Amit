import numpy as np
from astropy.io import fits as pyfits
from astropy.io import fits
import pandas
import sys
import fitsio
import glob
from sqlitedict import SqliteDict
import gnupg

gpg = gnupg.GPG()

def FITS_to_PGP(message):
    """
    Turns a string stored into an FITS comment back into a proper
    PGP message
    """
    s = "-----BEGIN PGP MESSAGE-----\n\n"
    s += message
    s += "\n-----END PGP MESSAGE-----\n"
    return s

# This is the function which you will use to decrypt the string
def decrypt_string(encrypted_string):
    #string = getpass.getpass(prompt="Enter passphrase to unlock the private key:\n")
    decrypted_data = gpg.decrypt(encrypted_string) #, passphrase=string)

    print 'ok: ', decrypted_data.ok
    print 'status: ', decrypted_data.status
    #print 'stderr: ', decrypted_data.stderr
    #print 'Decrypted string (additive m value for the catalog): ', decrypted_data.data
    return decrypted_data.data

def lens_select(lensargs):
    """
    lens_select(dictionary lensargs)

    Selects lens sample

    :Parameters:
    
    -    lensargs : Dictionary read from the config file which contains various
         cuts on the lens sample

    :Returns:
    
    -    ra : Right ascension in degrees
    -    dec : Declination in degrees
    -    zred : Redshift
    -    weight : Weights for the lenses

    :Adding new samples:
    
    Add new samples by copy pasting one of the lensargs block and then modifying them according to your will
    """
    if lensargs["type"] == "generic":
        hdulist = pyfits.open(lensargs["fname"])
        data = hdulist[1].data
        ra = data["ra"].astype("float64")
        dec = data["dec"].astype("float64")
        z = data["z"].astype("float64")
        wt = data["wt"].astype("float64")
         
        sys.stdout.write("Selecting %d samples \n" % (ra.size))

        return ra, dec, z, wt

    # Use generic text file with the columns: 
    # ra (deg), dec (deg), z, lens_wt
    # Any lines beginning with "#" are ignored
    if lensargs["type"] == "generic-text":
        data = pandas.read_csv(lensargs["fname"], delim_whitespace=1, header=None, names=(["ra", "dec", "z", "wt"]), comment="#")
        ra = data["ra"].values.astype("float64")
        dec = data["dec"].values.astype("float64")
        z = data["z"].values.astype("float64")
        wt = data["wt"].values.astype("float64")
         
        sys.stdout.write("Selecting %d samples \n" % (ra.size))

        return ra, dec, z, wt
    if lensargs["type"] == "gama":
        hdulist = pyfits.open(".//home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/gama/G3CFOFGroup.fits")
        data = hdulist[1].data
        ra = data["IterCenRA"].astype("float64") #ra of galaxies in the lense
        dec = data["IterCenDec"].astype("float64") #dec of galaxies in the lense
        z = data["IterCenZ"].astype("float64") #redshift of galaxies in the lense
        vdisp=data["VelDisp"].astype("float64") 
        idx = (z>lensargs["zmin"]) & (z<=lensargs["zmax"]) & (vdisp>lensargs["vdispmin"]) & (vdisp<=lensargs["vdispmax"])
        ra = ra[idx]
        dec = dec[idx]
        z = z[idx]
        vdisp = vdisp[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d samples \n" % (ra.size))

        return ra, dec, z, wt

    if lensargs["type"] == "mgii-absorbers-ran":
        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/INFO_140303_Ref_Mg_II_DR79_match.fits")
        data = hdulist[1].data
        ra = data["ra"].astype("float64")
        dec = data["dec"].astype("float64")
        zabs = data["zabs"].astype("float64")
         
        idx = (zabs>lensargs["zmin"]) & (zabs<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zabs = zabs[idx]

        wt = ra/ra

        sys.stdout.write("Selecting %d samples \n" % (ra.size))

        return ra, dec, zabs, wt

    if lensargs["type"] == "mgii-absorbers":
        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/INFO_140303_Mg_II_DR79_match.fits")
        data = hdulist[1].data
        ra = data["ra"].astype("float64")
        dec = data["dec"].astype("float64")
        zabs = data["zabs"].astype("float64")
         
        idx = (zabs>lensargs["zmin"]) & (zabs<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zabs = zabs[idx]

        wt = ra/ra

        sys.stdout.write("Selecting %d samples \n" % (ra.size))

        return ra, dec, zabs, wt

    if lensargs['type'] == "hsc-cmass-subsample-random":
        df = pandas.read_csv("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/HSC_CMASS/CMASS_K+E_v0.99.dat", delim_whitespace=1)
        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values
        wt = (df.weight_noz.values + df.weight_cp.values - 1.0)* df.weight_star.values
        iMag = df.Mag_Wake.values
        idx = (iMag>lensargs["iMagmin"]) & (iMag<=lensargs["iMagmax"]) & (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])

        # First compute the selection function as a function of redshift
        idxall = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        zarr = np.linspace(lensargs["zmin"], lensargs["zmax"], 20)
        selarr = np.zeros(zarr.size-1)
        for i in range(1, zarr.size):
            selarr[i-1] = np.sum(idx & (zred>zarr[i-1]) & (zred<=zarr[i]))*1.0/np.sum(idxall & (zred>zarr[i-1]) & (zred<=zarr[i]))
            print i, selarr[i-1], "\n================="
        zarr = zarr[:-1]/2. + zarr[1:]/2.

        from scipy.interpolate import UnivariateSpline
        selspl = UnivariateSpline(zarr, selarr, bbox=(lensargs["zmin"], lensargs["zmax"]), ext=1)

        with SqliteDict("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/HSC_CMASS/randoms/ran%03d.sqlite3" % (lensargs["rannum"])) as data:
            ra = data["ra"]
            dec = data["dec"]
            zred = data["z"].astype('float64')
            wt = data["weight"].astype("float64")

            idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
            ra = ra[idx]
            dec = dec[idx]
            zred = zred[idx]
            wt = wt[idx]

            np.random.seed(lensargs["rannum"])
            rannum = np.random.random(size=ra.size)
            idx = (rannum<selspl(zred))
            ra = ra[idx]
            dec = dec[idx]
            zred = zred[idx]
            wt = wt[idx]

        sys.stdout.write("Selecting %d samples \n" % (ra.size))
        #jackreg = getregions(ra, dec, lensargs["jackregfile"])
        return ra, dec, zred, wt #, jackreg

    if lensargs['type'] == "hsc-lowz-subsample":
        df = pandas.read_csv("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/HSC_LOWZ/LOWZ_K+E_v0.99.dat", delim_whitespace=1)
        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values
        wt = (df.weight_noz.values + df.weight_cp.values - 1.0)* df.weight_star.values
        iMag = df.Mag_Wake.values
        idx = (iMag>lensargs["iMagmin"]) & (iMag<=lensargs["iMagmax"]) & (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])

        if lensargs['selection_write'] == 1:
            idxall = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
            zarr = np.linspace(lensargs["zmin"], lensargs["zmax"], 20)
            selarr = np.zeros(zarr.size-1)
            for i in range(1, zarr.size):
                selarr[i-1] = np.sum(idx & (zred>zarr[i-1]) & (zred<=zarr[i]))*1.0/np.sum(idxall & (zred>zarr[i-1]) & (zred<=zarr[i]))
                print i, selarr[i-1], "\n================="
            zarr = zarr[:-1]/2. + zarr[1:]/2.
            dfsel = pandas.DataFrame(zarr, columns=(["zred"]))
            dfsel["sel"] = selarr
            dfsel.to_csv(lensargs["selection"], index=False, sep=" ")

        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
        iMag = iMag[idx]

        sys.stdout.write("Selecting %d samples \n" % (ra.size))
        #jackreg = getregions(ra, dec, lensargs["jackregfile"])
        return ra, dec, zred, wt #, jackreg

    if lensargs['type'] == "hsc-cmass-subsample":
        df = pandas.read_csv("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/HSC_CMASS/CMASS_K+E_v0.99.dat", delim_whitespace=1)
        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values
        wt = (df.weight_noz.values + df.weight_cp.values - 1.0)* df.weight_star.values
        iMag = df.Mag_Wake.values
        idx = (iMag>lensargs["iMagmin"]) & (iMag<=lensargs["iMagmax"]) & (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])

        if lensargs['selection_write'] == 1:
            idxall = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
            zarr = np.linspace(lensargs["zmin"], lensargs["zmax"], 20)
            selarr = np.zeros(zarr.size-1)
            for i in range(1, zarr.size):
                selarr[i-1] = np.sum(idx & (zred>zarr[i-1]) & (zred<=zarr[i]))*1.0/np.sum(idxall & (zred>zarr[i-1]) & (zred<=zarr[i]))
                print i, selarr[i-1], "\n================="
            zarr = zarr[:-1]/2. + zarr[1:]/2.
            dfsel = pandas.DataFrame(zarr, columns=(["zred"]))
            dfsel["sel"] = selarr
            dfsel.to_csv(lensargs["selection"], index=False, sep=" ")

        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
        iMag = iMag[idx]

        sys.stdout.write("Selecting %d samples \n" % (ra.size))
        #jackreg = getregions(ra, dec, lensargs["jackregfile"])
        return ra, dec, zred, wt #, jackreg

    if lensargs['type'] == "SDSS-QSOs-John":
        fname = "/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/SDSS_QSOs_zlt1.asc"
        dfreal = pandas.read_csv(fname, delim_whitespace=1, names=(["ra", "dec", "z"]), usecols=([1, 2, 7]))
        ra = dfreal["ra"].values
        dec = dfreal["dec"].values
        zred = dfreal["z"].values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt

    if lensargs['type'] == "Masato-mock-forreq-random-full":
        rannum = lensargs['rannum']
        rotation = lensargs['rotation']
        realization = lensargs['realization']
        files = glob.glob("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Mocks_HOD_full/random_*")
        np.random.seed(rannum)

        for i, fname in enumerate(files):
            tractreg = fname.split("_")[-1]
            realfile = "/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Mocks_HOD_full/r%03d_rotmat%d_%s" % (realization, rotation, tractreg)
            dfreal = pandas.read_csv(realfile, delim_whitespace=1, skiprows=1, names=(["ra", "dec", "z"]), usecols=([1, 2, 3]))
            Nreal = dfreal.ra.values.size

            rowstart = rannum * Nreal
            if i==0:
                df = pandas.read_csv(fname, delim_whitespace=1, skiprows=1+rowstart, names=(["ra", "dec", "z"]), usecols=([0, 1, 2]), nrows=Nreal)
            else:
                dfp = pandas.read_csv(fname, delim_whitespace=1, skiprows=1+rowstart, names=(["ra", "dec", "z"]), usecols=([0, 1, 2]), nrows=Nreal)
                df = df.append(dfp, ignore_index=1)

        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt


    if lensargs['type'] == "Masato-mock-forreq-random":
        rannum = lensargs['rannum']
        rotation = lensargs['rotation']
        realization = lensargs['realization']
        files = glob.glob("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Mocks_HOD/random_*")
        np.random.seed(rannum)

        for i, fname in enumerate(files):
            tractreg = fname.split("_")[-1]
            realfile = "/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Mocks_HOD/r%03d_rotmat%d_%s" % (realization, rotation, tractreg)
            dfreal = pandas.read_csv(realfile, delim_whitespace=1, skiprows=1, names=(["ra", "dec", "z"]), usecols=([1, 2, 3]))
            Nreal = dfreal.ra.values.size

            rowstart = rannum * Nreal
            if i==0:
                df = pandas.read_csv(fname, delim_whitespace=1, skiprows=1+rowstart, names=(["ra", "dec", "z"]), usecols=([0, 1, 2]), nrows=Nreal)
            else:
                dfp = pandas.read_csv(fname, delim_whitespace=1, skiprows=1+rowstart, names=(["ra", "dec", "z"]), usecols=([0, 1, 2]), nrows=Nreal)
                df = df.append(dfp, ignore_index=1)

        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt

    if lensargs['type'] == "Masato-mock-forreq-full":
        realization = lensargs['realization']
        rotation = lensargs['rotation']
        files = glob.glob("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Mocks_HOD_full/r%03d_rotmat%d_*" % (realization, rotation))
        for i, fname in enumerate(files):
            if i==0:
                df = pandas.read_csv(fname, delim_whitespace=1, skiprows=1, names=(["id", "ra", "dec", "z"]), usecols=([0, 1, 2, 3]))
            else:
                df = df.append(pandas.read_csv(fname, delim_whitespace=1, skiprows=1, names=(["id", "ra", "dec", "z"]), usecols=([0, 1, 2, 3])), ignore_index=1)
        
        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt


    if lensargs['type'] == "Masato-mock-forreq":
        realization = lensargs['realization']
        rotation = lensargs['rotation']
        files = glob.glob("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Mocks_HOD/r%03d_rotmat%d_*" % (realization, rotation))
        for i, fname in enumerate(files):
            if i==0:
                df = pandas.read_csv(fname, delim_whitespace=1, skiprows=1, names=(["id", "ra", "dec", "z"]), usecols=([0, 1, 2, 3]))
            else:
                df = df.append(pandas.read_csv(fname, delim_whitespace=1, skiprows=1, names=(["id", "ra", "dec", "z"]), usecols=([0, 1, 2, 3])), ignore_index=1)
        
        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt

    if lensargs['type'] == "Masato-mock":
        realization = lensargs['realization']
        rotation = lensargs['rotation']
        sample = lensargs['sample']
        zmin = lensargs["zmin"]
        zmax = lensargs["zmax"]
        files = glob.glob("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/hod_gal/%s_%.2f_%.2f/newID/r%03d_rotmat%d_*" % (sample, zmin, zmax, realization, rotation))
        for i, fname in enumerate(files):
            if i==0:
                df = pandas.read_csv(fname, delim_whitespace=1, skiprows=1, names=(["id", "ra", "dec", "z"]), usecols=([0, 1, 2, 3]))
            else:
                df = df.append(pandas.read_csv(fname, delim_whitespace=1, skiprows=1, names=(["id", "ra", "dec", "z"]), usecols=([0, 1, 2, 3])), ignore_index=1)
        
        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt

    if lensargs['type'] == "dr12-lowz":
        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/DR12/galaxy_DR12v5_LOWZ_North.fits.gz")
        data = hdulist[1].data
        ra_n = data["ra"]
        dec_n = data["dec"]
        zred_n = data["z"].astype('float64')
        wt_n = data["weight_systot"].astype("float64")
        idx = (zred_n>lensargs["zmin"]) & (zred_n<=lensargs["zmax"])
        ra_n = ra_n[idx]
        dec_n = dec_n[idx]
        zred_n = zred_n[idx]
        wt_n = wt_n[idx]

        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/DR12/galaxy_DR12v5_LOWZ_South.fits.gz")
        data = hdulist[1].data
        ra_s = data["ra"]
        dec_s = data["dec"]
        zred_s = data["z"].astype('float64')
        wt_s = data["weight_systot"].astype("float64")
        idx = (zred_s>lensargs["zmin"]) & (zred_s<=lensargs["zmax"])
        ra_s = ra_s[idx]
        dec_s = dec_s[idx]
        zred_s = zred_s[idx]
        wt_s = wt_s[idx]

        ra = np.append(ra_n, ra_s)
        dec = np.append(dec_n, dec_s)
        zred = np.append(zred_n, zred_s)
        wt = np.append(wt_n, wt_s)
        
        sys.stdout.write("Selecting %d lenses \n" % (ra.size))

        if "pofzoutput" in lensargs:
            n, bins = np.histogram(zred, weights=wt, bins=np.linspace(np.min(zred), np.max(zred), 20), normed=1)
            np.savetxt(lensargs["pofzoutput"], np.transpose([bins[1:]/2+bins[:-1]/2, n]))
            exit(11)

        return ra, dec, zred, wt

    if lensargs['type'] == "cfht-cmass-reproduce":
        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/HSC_CMASS/cmass-dr11v1-A-Anderson.dat_wcollind_mstel_xdist.fits")
        data = hdulist[1].data
        ra = data["ra"]
        dec = data["dec"]
        zred = data["z"].astype('float64')
        wt_star = data["weight_star"].astype("float64")
        wt_noz = data["weight_noz"].astype("float64")
        wt_cp = data["weight_cp"].astype("float64")
        wt = wt_star*(wt_noz+wt_cp-1)
        mstar = data["logmass"]
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (mstar>lensargs["mstarmin"]) & (mstar<=lensargs["mstarmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
        
        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt

    if lensargs['type'] == "dr12-cmass":
        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/DR12/galaxy_DR12v5_CMASS_North.fits.gz")
        data = hdulist[1].data
        ra_n = data["ra"]
        dec_n = data["dec"]
        zred_n = data["z"].astype('float64')
        wt_n = data["weight_systot"].astype("float64")
        idx = (zred_n>lensargs["zmin"]) & (zred_n<=lensargs["zmax"])
        ra_n = ra_n[idx]
        dec_n = dec_n[idx]
        zred_n = zred_n[idx]
        wt_n = wt_n[idx]

        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/DR12/galaxy_DR12v5_CMASS_South.fits.gz")
        data = hdulist[1].data
        ra_s = data["ra"]
        dec_s = data["dec"]
        zred_s = data["z"].astype('float64')
        wt_s = data["weight_systot"].astype("float64")
        idx = (zred_s>lensargs["zmin"]) & (zred_s<=lensargs["zmax"])
        ra_s = ra_s[idx]
        dec_s = dec_s[idx]
        zred_s = zred_s[idx]
        wt_s = wt_s[idx]

        ra = np.append(ra_n, ra_s)
        dec = np.append(dec_n, dec_s)
        zred = np.append(zred_n, zred_s)
        wt = np.append(wt_n, wt_s)
        
        sys.stdout.write("Selecting %d lenses \n" % (ra.size))

        if "pofzoutput" in lensargs:
            n, bins = np.histogram(zred, weights=wt, bins=np.linspace(np.min(zred), np.max(zred), 20), normed=1)
            np.savetxt(lensargs["pofzoutput"], np.transpose([bins[1:]/2+bins[:-1]/2, n]))
            exit(11)

        return ra, dec, zred, wt

    if lensargs['type'] == "voids":

        df = pandas.read_csv("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Voids/Void_combined_catalog.dat", delim_whitespace=1)
        ra = df.ra.values
        dec = df.dec.values
        zred = df.zred.values
        wt = ra/ra
        
        sys.stdout.write("Selecting %d voids \n" % (ra.size))
        #return ra[idx], dec[idx], zred[idx], wt[idx]
        return ra, dec, zred, wt

    if lensargs['type'] == "alexie-test":
        ra, dec, zred, wt = np.loadtxt("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/alexie-test/bossdr12_GAMA15_zbin1.txt", unpack=1)
        wt = ra/ra
        
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        #return ra[idx], dec[idx], zred[idx], wt[idx]
        return ra, dec, zred, wt

    if lensargs['type'] == "redmapper-smallR":
        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/redmapper-smallR/dr8_run_redmapper_v5.10-0.75_lgt5_prunedcatalog.fit")
        data = hdulist[1].data
        ra = data["ra"]
        dec = data["dec"]
        zred = data["z_lambda"].astype('float64')
        lamda = data["lambda_chisq"]
        wt = ra/ra
        
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lamda>lensargs["lammin"]) & (lamda<=lensargs["lammax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        #return ra[idx], dec[idx], zred[idx], wt[idx]
        return ra, dec, zred, wt

    if lensargs['type'] == "redmapper-largeR":
        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/redmapper-largeR/dr8_run_redmapper_v5.10-1.25_lgt5_prunedcatalog.fit")
        data = hdulist[1].data
        ra = data["ra"]
        dec = data["dec"]
        zred = data["z_lambda"].astype('float64')
        lamda = data["lambda_chisq"]
        wt = ra/ra
        
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lamda>lensargs["lammin"]) & (lamda<=lensargs["lammax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        #return ra[idx], dec[idx], zred[idx], wt[idx]
        return ra, dec, zred, wt

    if lensargs['type'] == "redmapper":
        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/redmapper/redmapper_dr8_public_v5.10_catalog.fits")
        data = hdulist[1].data
        ra = data["ra"]
        dec = data["dec"]
        zred = data["z_lambda"].astype('float64')
        lamda = data["lambda"]
        wt = ra/ra
        
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lamda>lensargs["lammin"]) & (lamda<=lensargs["lammax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        #return ra[idx], dec[idx], zred[idx], wt[idx]
        return ra, dec, zred, wt

    if lensargs['type'] == "arindam_lenses":
        fp=pandas.read_csv("arindam/arindam_lenses.csv",delimiter=",",usecols=[0,1,2],header=None,skiprows=1)
        ra = fp[0].astype('float64')
        dec = fp[1].astype('float64')
        zred = fp[2].astype('float64')
        wt = ra/ra
        idx=fp[2] > 0
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        return ra, dec, zred, wt



    if lensargs['type'] == "redmapper_satellite":
        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/redmapper/redmapper_dr8_public_v6.3_members.fits")
        hdulist2 = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/redmapper/redmapper_dr8_public_v6.3_catalog.fits")
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
        data=data[tidx]



        zlamda=data2["z_lambda"]
        cid=data2["id"]
        ldict={}
        for ii in range(cid.size):
            ldict[cid[ii]]=zlamda[ii]
        mid=data['id']
        ra = data["ra"]
        dec = data["dec"]        
        lamda = data2["lambda"]
        wt = ra/ra
        r = data["r"]
        zred = ra*0.0
        
        for jj in range(mid.size):
            zred[jj]=ldict[mid[jj]].astype('float64')        
        
        rr=r*(1.+zred)   #to convert from physical to comoving distances
    
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (rr>lensargs["rmin"]) & (rr<=lensargs["rmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        #lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        #return ra[idx], dec[idx], zred[idx], wt[idx]
        return ra, dec, zred, wt

    
    if lensargs['type'] == "redmapper_satellite_stellar_counterpart":
        #hdulist = pyfits.open("pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/redmapper/redmapper_dr8_public_v6.3_members.fits")
        hdulist=pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/redmapper/redmapper_hsc_sdss_stellar_counterpart.fits")        
        #hdulist2 = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/redmapper/redmapper_dr8_public_v6.3_catalog.fits")
        data = hdulist[1].data
        #data2 = hdulist2[1].data
        #zlamda=data2["z_lambda"]
        #cid=data2["id"]
        #ldict={}
        #for ii in range(cid.size):
        #    ldict[cid[ii]]=zlamda[ii]
        #mid=data['id']
        ra = data["sdssdata_ra"].astype("float64")
        dec = data["sdssdata_dec"].astype("float64")        
        #lamda = data2["lambda"].astype("float64")
        wt = ra/ra
        r = data["redmapper_rsat"]
        zred = data["redmapper_zstellar"].astype("float64")
        
        #for jj in range(mid.size):
        #    zred[jj]=ldict[mid[jj]].astype('float64')        
            
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (r>lensargs["rmin"]) & (r<=lensargs["rmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        #lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        #return ra[idx], dec[idx], zred[idx], wt[idx]
        return ra, dec, zred, wt


    if lensargs['type'] == "redmapper_random":
        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/redmapper/redmapper_dr8_public_v6.3_randoms.fits")
        data = hdulist[1].data
        ra = data["ra"].astype('float64')
        dec = data["dec"].astype('float64')
        zred = data["z"].astype('float64')
        lamda = data["lambda"].astype('float64')
        wt = data["weight"].astype('float64')
        
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        #return ra[idx], dec[idx], zred[idx], wt[idx]
        return ra, dec, zred, wt




 	
    if lensargs['type'] == "redmapper-random":
        ra, dec, zred, lamda, wt = np.loadtxt("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/redmapper/redmapper_public_v5.10_randoms_%05d.dat" % (lensargs["rannum"]), unpack=1) 
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lamda>lensargs["lammin"]) & (lamda<=lensargs["lammax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d randoms \n" % (np.sum(idx)))
        return ra[idx], dec[idx], zred[idx], wt[idx]




    if lensargs['type'] == "redmapper_satellite_in_hsc_field":
        hdulist = pyfits.open("hsc_redmapper_overlap/redmapper_members_hsc_overlap_.fits")
        data = hdulist[1].data
        
        ra = data['redmapper_Ra'].astype('float64')
        dec = data['redmapper_Dec'].astype('float64')
        wt = ra/ra
        r = data['dist_from_cluster_center'].astype('float64')
        zred = data['redmapper_redshift'].astype('float64')
        
        
        rr=r*(1.+zred)   #to convert from physical to comoving distances
    
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (rr>lensargs["rmin"]) & (rr<=lensargs["rmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        #lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        #return ra[idx], dec[idx], zred[idx], wt[idx]
        return ra, dec, zred, wt




    if lensargs['type'] == "sdss-vagc-lbg":
        ra, dec, zred, mstel, mag = np.loadtxt("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/LBG/Isolated_v0.dat", unpack=1)
        mstel = np.log10(mstel)
        wt = ra/ra

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (mstel>lensargs["xmstelmin"]) & (mstel<=lensargs["xmstelmax"])
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        return ra[idx], dec[idx], zred[idx], wt[idx]

    if lensargs['type'] == "lowz_rcsrestrict":
        ra, dec, zred, lam, wt = np.loadtxt("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/rcsrestrict_dr12_v5_lowz_All_catalog.dat", unpack=1)

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        return ra[idx], dec[idx], zred[idx], wt[idx]

    if lensargs['type'] == "camira":
        ra, dec, zred, lam = np.loadtxt("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/camira/%s/%s/camira_catalog.dat" % (lensargs['hsc-release'], lensargs['version']), unpack=1, usecols=(0, 1, 2, 3))
        wt = lam/lam

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lam>lensargs["lammin"]) & (lam<=lensargs["lammax"])
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))

        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]

        if "pofzoutput" in lensargs:
            n, bins = np.histogram(zred, weights=wt, bins=np.linspace(np.min(zred), np.max(zred), 20), normed=1)
            np.savetxt(lensargs["pofzoutput"], np.transpose([bins[1:]/2+bins[:-1]/2, n]))
            exit(11)
        return ra, dec, zred, wt


def source_select(sourceargs, chunksize):
    """
    source_select(dictionary lensargs)

    Selects source sample

    :Parameters:
    
    -    sourceargs : Dictionary read from the config file which contains various
         cuts on the source sample
    -    chunksize : The chunksize of sources to be read (only applicable for some ascii readers)

    :Returns:     return datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag
    
    -    datagal : array consisting of ra, dec, e1, e2, weight, erms, m, c1, c2, a1, a2
    -    sourceargs : Dictionary read from the config file which contains various
         cuts on the source sample
    -    Ngal : Number of galaxies
    -    status : status of the read
    -    pofz : Record arrays read from the P(z) file for `Ngal` galaxies
    -    datagalflag: Flag array to indicate whether to use the galaxy or not
    -    pofzflag: Flag array to indicate whether to use the photoz of the galaxy or not

    """
    if sourceargs['type'] == "hsc-wide-s16a_v2.0" and sourceargs['filetype'] == "fits" and not sourceargs["fullpofz"]:

        # This variable iterates over the different fields in the catalog
        itern = sourceargs['iter']
        if itern == 0:
            if not sourceargs["blinded"]:
                sourceargs["df"] = pandas.read_csv("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/Mapping_from_Rachel_catalogs_to_Hironao_catalogs.dat.calibrated", delim_whitespace=1)
            else:
                sourceargs["df"] = pandas.read_csv("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/Mapping_from_Rachel_catalogs_to_Hironao_catalogs.dat", delim_whitespace=1)
            if "usefield" in sourceargs:
                restidx = sourceargs["df"].field.values==sourceargs["usefield"]
                sourceargs["df"] = sourceargs["df"][restidx]
                print "Using only field", sourceargs["usefield"]

        datagal = 0
        Ngal = 0
        pofz = 0
        status = not (itern<(sourceargs['df'].field.values.size))
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0., 0.

        field = sourceargs["df"].field.values[itern]
        tract = sourceargs["df"].tract.values[itern]
        begin = sourceargs["df"].begin.values[itern]
        end = sourceargs["df"].end.values[itern]
        
        # Initialize list of tracts within a given field
        if begin == 0:
            # Define the shape catalog file from Rachel to read from for this field
            if sourceargs["blinded"]:
                # Decrypt the Delta m value, if sourceargs["dm"] not set already
                fname = "/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Blinded_S16A_v2.0/%s_blinded_m_%s_%d.fits" % (field, sourceargs["username"], sourceargs["blindnumber"])
                fname_nom = "/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Blinded_S16A_v2.0/%s_no_m.fits" % (field)
                if "dm" not in sourceargs:
                    try:
                        hdulist = pyfits.open(fname, memmap=True)
                        msg = hdulist[0].header["DM1"]
                        sourceargs['dm'] = float(decrypt_string(FITS_to_PGP(msg)))
                    except:
                        print "GPG decryption failed, check your gpg-agent"
                        exit(12)
                # The commented line is memory intensive
                #sourceargs["fits"] = fitsio.FITS("%s[1][col *, ishape_hsm_regauss_derived_shear_bias_m = ishape_hsm_regauss_derived_shear_bias_m - %f ]" % (fname, sourceargs['dm']))[1]
                # Instead read the file as it is and then subtract off later
                sourceargs["fits_mblind"] = fitsio.FITS("%s" % (fname))[1]
                sourceargs["fits"] = fitsio.FITS("%s" % (fname_nom))[1]
            else:
                fname = "/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Calibrated_S16A_v2.0/%s_calibrated.fits" % (field)
                sourceargs["fits"] = fitsio.FITS("%s" % (fname))[1]
                sourceargs["fits_mblind"] = fitsio.FITS("%s" % (fname))[1]
                sourceargs["dm"] = 0.0

        # Read in the photoz file for the current field, tract
        try:
            hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/%s_tracts/%d_pz.fits" % (field, tract))
            data = hdulist[1].data
            pofz = data["%s_photoz%s" % (sourceargs['pofz_type'], sourceargs['photoz_estimate']) ]
            pofzflag = (data["%s_photoz_risk_best" % (sourceargs['pofz_type']) ]<sourceargs['photoz_risk_best_cut'])
        except:
            pofz = 0
            pofzflag = np.zeros(end-begin, dtype=bool)

        # Read in the shape values for the current field, tract, the row filtering specification indexes rows from 1 - Nrecords
        wheremask = sourceargs['fits'].where("#row>%d && #row<=%d" % (begin, end))

        #datagal = sourceargs['fits'][1]["ira", "idec", "ishape_hsm_regauss_e1", "ishape_hsm_regauss_e2", "ishape_hsm_regauss_derived_shear_bias_c1", "ishape_hsm_regauss_derived_shear_bias_c2", "ishape_hsm_regauss_derived_shear_bias_m", "ishape_hsm_regauss_derived_shape_weight", "ishape_hsm_regauss_derived_rms_e"][wheremask]
        #datagal = sourceargs['fits']["ira", "idec", "ishape_hsm_regauss_e1", "ishape_hsm_regauss_e2", "ishape_hsm_regauss_derived_shape_weight", "ishape_hsm_regauss_derived_rms_e", "ishape_hsm_regauss_derived_shear_bias_m", "ishape_hsm_regauss_derived_shear_bias_c1", "ishape_hsm_regauss_derived_shear_bias_c2"][wheremask]
        ira = sourceargs['fits']["ira"][wheremask]
        idec = sourceargs['fits']["idec"][wheremask]
        ishape_hsm_regauss_e1 = sourceargs['fits']["ishape_hsm_regauss_e1"][wheremask]
        ishape_hsm_regauss_e2 = sourceargs['fits']["ishape_hsm_regauss_e2"][wheremask]
        ishape_hsm_regauss_derived_shape_weight = sourceargs['fits']["ishape_hsm_regauss_derived_shape_weight"][wheremask]
        ishape_hsm_regauss_derived_rms_e = sourceargs['fits']["ishape_hsm_regauss_derived_rms_e"][wheremask]
        ishape_hsm_regauss_derived_shear_bias_m = sourceargs['fits_mblind']["ishape_hsm_regauss_derived_shear_bias_m"][wheremask] - sourceargs["dm"]
        ishape_hsm_regauss_derived_shear_bias_c1 = sourceargs['fits']["ishape_hsm_regauss_derived_shear_bias_c1"][wheremask]
        ishape_hsm_regauss_derived_shear_bias_c2 = sourceargs['fits']["ishape_hsm_regauss_derived_shear_bias_c2"][wheremask]
        ishape_hsm_regauss_resolution = sourceargs['fits']["ishape_hsm_regauss_resolution"][wheremask]
        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_derived_shape_weight, ishape_hsm_regauss_derived_rms_e, ishape_hsm_regauss_derived_shear_bias_m, ishape_hsm_regauss_derived_shear_bias_c1, ishape_hsm_regauss_derived_shear_bias_c2, ishape_hsm_regauss_resolution]).T
        datagalflag = sourceargs['fits']["weak_lensing_flag"][wheremask]

        sourceargs['iter'] = sourceargs['iter'] + 1

        Ngal = pofzflag.size
        status = 0
        
        return datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag

    if sourceargs['type'] == "hsc-wide-s16a_v2.0" and sourceargs['filetype'] == "fits" and  sourceargs["fullpofz"]:

        # This variable iterates over the different fields in the catalog
        itern = sourceargs['iter']
        if itern == 0:
            if not sourceargs["blinded"]:
                sourceargs["df"] = pandas.read_csv("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/Mapping_from_Rachel_catalogs_to_Hironao_catalogs.dat.calibrated", delim_whitespace=1)
            else:
                sourceargs["df"] = pandas.read_csv("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/Mapping_from_Rachel_catalogs_to_Hironao_catalogs.dat", delim_whitespace=1)
            if "usefield" in sourceargs:
                restidx = sourceargs["df"].field.values==sourceargs["usefield"]
                sourceargs["df"] = sourceargs["df"][restidx]
                print "Using only field", sourceargs["usefield"]

        datagal = 0
        Ngal = 0
        pofz = 0
        status = not (itern<(sourceargs['df'].field.values.size))
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0., 0.

        field = sourceargs["df"].field.values[itern]
        tract = sourceargs["df"].tract.values[itern]
        begin = sourceargs["df"].begin.values[itern]
        end = sourceargs["df"].end.values[itern]
        
        # Initialize list of tracts within a given field
        if begin == 0:
            # Define the shape catalog file from Rachel to read from for this field
            if sourceargs["blinded"]:
                # Decrypt the Delta m value, if sourceargs["dm"] not set already
                fname = "/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Blinded_S16A_v2.0/%s_blinded_m_%s_%d.fits" % (field, sourceargs["username"], sourceargs["blindnumber"])
                fname_nom = "/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Blinded_S16A_v2.0/%s_no_m.fits" % (field)
                if "dm" not in sourceargs:
                    try:
                        hdulist = pyfits.open(fname, memmap=True)
                        msg = hdulist[0].header["DM1"]
                        sourceargs['dm'] = float(decrypt_string(FITS_to_PGP(msg)))
                    except:
                        print "GPG decryption failed, check your gpg-agent"
                        exit(12)
                # The commented line is memory intensive
                #sourceargs["fits"] = fitsio.FITS("%s[1][col *, ishape_hsm_regauss_derived_shear_bias_m = ishape_hsm_regauss_derived_shear_bias_m - %f ]" % (fname, sourceargs['dm']))[1]
                # Instead read the file as it is and then subtract off later
                sourceargs["fits_mblind"] = fitsio.FITS("%s" % (fname))[1]
                sourceargs["fits"] = fitsio.FITS("%s" % (fname_nom))[1]
            else:
                fname = "/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Calibrated_S16A_v2.0/%s_calibrated.fits" % (field)
                sourceargs["fits"] = fitsio.FITS("%s" % (fname))[1]
                sourceargs["fits_mblind"] = fitsio.FITS("%s" % (fname))[1]
                sourceargs["dm"] = 0.0

        # Read in the photoz file for the current field, tract
        try:
            hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/%s_tracts/%d_pz.fits" % (field, tract))
            data = hdulist[1].data
            pofz = data["%s_photoz_best" % (sourceargs['pofz_type']) ]
            pofzflag = (data["%s_photoz_risk_best" % (sourceargs['pofz_type']) ]<sourceargs['photoz_risk_best_cut'])
        except:
            pofz = 0
            pofzflag = np.zeros(end-begin, dtype=bool)

        # Read in the shape values for the current field, tract
        wheremask = sourceargs['fits'].where("#row>%d && #row<=%d" % (begin, end))

        #datagal = sourceargs['fits']["ira", "idec", "ishape_hsm_regauss_e1", "ishape_hsm_regauss_e2", "ishape_hsm_regauss_derived_shear_bias_c1", "ishape_hsm_regauss_derived_shear_bias_c2", "ishape_hsm_regauss_derived_shape_weight", "ishape_hsm_regauss_derived_rms_e"][wheremask]
        #datagalflag = sourceargs['fits'][1]["weak_lensing_flag"][wheremask]

        ira = sourceargs['fits']["ira"][wheremask]
        idec = sourceargs['fits']["idec"][wheremask]
        ishape_hsm_regauss_e1 = sourceargs['fits']["ishape_hsm_regauss_e1"][wheremask]
        ishape_hsm_regauss_e2 = sourceargs['fits']["ishape_hsm_regauss_e2"][wheremask]
        ishape_hsm_regauss_derived_shape_weight = sourceargs['fits']["ishape_hsm_regauss_derived_shape_weight"][wheremask]
        ishape_hsm_regauss_derived_rms_e = sourceargs['fits']["ishape_hsm_regauss_derived_rms_e"][wheremask]
        ishape_hsm_regauss_derived_shear_bias_m = sourceargs['fits_mblind']["ishape_hsm_regauss_derived_shear_bias_m"][wheremask] - sourceargs["dm"]
        ishape_hsm_regauss_derived_shear_bias_c1 = sourceargs['fits']["ishape_hsm_regauss_derived_shear_bias_c1"][wheremask]
        ishape_hsm_regauss_derived_shear_bias_c2 = sourceargs['fits']["ishape_hsm_regauss_derived_shear_bias_c2"][wheremask]
        ishape_hsm_regauss_resolution = sourceargs['fits']["ishape_hsm_regauss_resolution"][wheremask]
        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_derived_shape_weight, ishape_hsm_regauss_derived_rms_e, ishape_hsm_regauss_derived_shear_bias_m, ishape_hsm_regauss_derived_shear_bias_c1, ishape_hsm_regauss_derived_shear_bias_c2, ishape_hsm_regauss_resolution]).T
        datagalflag = sourceargs['fits']["weak_lensing_flag"][wheremask]

        try:
            sourceargs['fitspofz'] = fitsio.FITS("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/%s_tracts/%d_pz_pdf_%s.fits" % (field, tract, sourceargs['pofz_type']))[1]

            pofz = sourceargs['fitspofz']['P(z)'][:]
            pofz = pofz.reshape((pofz.shape[0], -1,))
        except:
            pofz = 0
            pofzflag = np.zeros(end-begin, dtype=bool)

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = pofzflag.size
        status = 0

        return datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag


    # S16A_v2.0 kept here for notes when required. This was before Rachel computed the calibration factors
    if sourceargs['type'] == "hsc-wide-s16a_v2.0_old" and sourceargs['filetype'] == "fits" and not sourceargs["fullpofz"]:

        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits_list'] = glob.glob("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/*_tracts/*_pz.fits")
            
        if itern == len(sourceargs['fits_list']):
            status = 1

        datagal = 0
        status = not (itern<len(sourceargs['fits_list']))
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0., 0.

        hdulist = pyfits.open(sourceargs['fits_list'][itern].split("_pz")[0] + ".fits")
        data = hdulist[1].data
        ira = data["ira"]
        idec = data["idec"]

        ishape_hsm_regauss_e1 = data["ishape_hsm_regauss_e1"]
        ishape_hsm_regauss_e2 = data["ishape_hsm_regauss_e2"]
        ishape_hsm_regauss_derived_bias_c1 = data["ishape_hsm_regauss_derived_bias_c1"]
        ishape_hsm_regauss_derived_bias_c2 = data["ishape_hsm_regauss_derived_bias_c2"]
        ishape_hsm_regauss_derived_bias_m = data["ishape_hsm_regauss_derived_bias_m"]
        ishape_hsm_regauss_derived_weight = data["ishape_hsm_regauss_derived_weight"] 
        ishape_hsm_regauss_derived_rms_e = data["ishape_hsm_regauss_derived_rms_e"] 
        zeroarr = ishape_hsm_regauss_derived_bias_c1 * 0.0
        datagalflag = data["weak_lensing_flag"]

        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_derived_weight, ishape_hsm_regauss_derived_rms_e, ishape_hsm_regauss_derived_bias_m, ishape_hsm_regauss_derived_bias_c1, ishape_hsm_regauss_derived_bias_c2]).T

        hdulist = pyfits.open(sourceargs['fits_list'][itern])
        data = hdulist[1].data
        pofz = data["%s_photoz_best" % (sourceargs['pofz_type']) ]
        pofzflag = (data["%s_photoz_risk_best" % (sourceargs['pofz_type']) ]<sourceargs['photoz_risk_best_cut'])

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = pofz.size
        status = 0
        return datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag

    '''
    if sourceargs['type'] == "hsc-wide-s16a_v2.0" and sourceargs['filetype'] == "fits" and  sourceargs["fullpofz"]:
        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits_list'] = glob.glob("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/*_tracts/*_%s.fits" % (sourceargs['pofz_type']))
            
        if itern == len(sourceargs['fits_list']):
            status = 1

        datagal = 0
        status = not (itern<len(sourceargs['fits_list']))
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0., 0.

        tractnum = (sourceargs['fits_list'][itern].split("/")[-1]).split("_pz")[0]
        tractreg = (sourceargs['fits_list'][itern].split("/")[2]).split("_tracts")[0]

        hdulist = pyfits.open(sourceargs['fits_list'][itern].split("_pz_pdf_%s" % (sourceargs['pofz_type']))[0] + ".fits")
        data = hdulist[1].data
        ira = data["ira"]
        idec = data["idec"]

        ishape_hsm_regauss_e1 = data["ishape_hsm_regauss_e1"]
        ishape_hsm_regauss_e2 = data["ishape_hsm_regauss_e2"]
        ishape_hsm_regauss_derived_bias_c1 = data["ishape_hsm_regauss_derived_bias_c1"]
        ishape_hsm_regauss_derived_bias_c2 = data["ishape_hsm_regauss_derived_bias_c2"]
        ishape_hsm_regauss_derived_bias_m = data["ishape_hsm_regauss_derived_bias_m"]
        ishape_hsm_regauss_derived_weight = data["ishape_hsm_regauss_derived_weight"] 
        ishape_hsm_regauss_derived_rms_e = data["ishape_hsm_regauss_derived_rms_e"] 
        zeroarr = ishape_hsm_regauss_derived_bias_c1 * 0.0
        datagalflag = data["weak_lensing_flag"]

        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_derived_weight, ishape_hsm_regauss_derived_rms_e, ishape_hsm_regauss_derived_bias_m, ishape_hsm_regauss_derived_bias_c1, ishape_hsm_regauss_derived_bias_c2, zeroarr, zeroarr]).T

        hdulist = pyfits.open(sourceargs['fits_list'][itern].split("_pdf_%s" % sourceargs['pofz_type'])[0] + ".fits")
        data = hdulist[1].data
        if sourceargs['pofz_type'] == "frankenz":
            pofzflag = (datagalflag == True)
        elif sourceargs['pofz_type'] == "mizuki":
            pofzflag = (datagalflag == True)
        elif sourceargs['pofz_type'] == "mlz":
            pofzflag = (datagalflag == True)
        elif sourceargs['pofz_type'] == "ephor":
            pofzflag = (datagalflag == True)

        sourceargs['fitspofz'] = fitsio.FITS(sourceargs['fits_list'][itern])

        pofz = sourceargs['fitspofz'][1]['P(z)'][:]
        pofz = pofz.reshape((pofz.shape[0], -1,))

        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/%s_tracts/%s_pz.fits" % (tractreg, tractnum) )
        data = hdulist[1].data
        pofzflag = (data["%s_photoz_risk_best" % (sourceargs['pofz_type']) ]<0.5)

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]
        status = 0

    if sourceargs['type'] == "mockhsc-wide-s16a_v2.0" and sourceargs['filetype'] == "fits" and  sourceargs["fullpofz"]:
        rotation = sourceargs['rotation']
        realization = sourceargs['realization']
        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits_list'] = glob.glob("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/*_tracts/*_%s.fits" % (sourceargs['pofz_type']))
            #sourceargs['fits_list_2'] = glob.glob("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/*_tracts/*_pz.fits")

        while itern<len(sourceargs['fits_list']):
            tractnum = (sourceargs['fits_list'][itern].split("/")[-1]).split("_pz")[0]
            tractreg = (sourceargs['fits_list'][itern].split("/")[2]).split("_tracts")[0]
            print tractnum, tractreg
            try:
                mockname = "/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/Mocks_shape/%s/r%03d/rotmat%s/mock_%s.fits" % (tractreg, realization, rotation, tractnum)
                print mockname
                hdulistmock = pyfits.open(mockname)
                datamock = hdulistmock[1].data
                break
            except:
                itern = itern + 1
                sourceargs['iter'] = itern
        
        if itern == len(sourceargs['fits_list']):
            status = 1

        datagal = 0
        status = not (itern<len(sourceargs['fits_list']))
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0., 0.

        print itern, sourceargs['fits_list'][itern]
        shapefile = sourceargs['fits_list'][itern].split("_pz_pdf_%s" % (sourceargs['pofz_type']))[0] + ".fits"
        print shapefile
        hdulist = pyfits.open(shapefile)
        data = hdulist[1].data
        ira = data["ira"]
        idec = data["idec"]

        ishape_hsm_regauss_e1 = datamock["e1_mock"].astype("float64")
        ishape_hsm_regauss_e2 = datamock["e2_mock"].astype("float64")

        ishape_hsm_regauss_derived_bias_c1 = data["ishape_hsm_regauss_derived_bias_c1"] * 0.0
        ishape_hsm_regauss_derived_bias_c2 = data["ishape_hsm_regauss_derived_bias_c2"] * 0.0
        ishape_hsm_regauss_derived_bias_m = data["ishape_hsm_regauss_derived_bias_m"] * 0.0
        ishape_hsm_regauss_derived_weight = data["ishape_hsm_regauss_derived_weight"] 
        ishape_hsm_regauss_derived_rms_e = data["ishape_hsm_regauss_derived_rms_e"] 
        zeroarr = ishape_hsm_regauss_derived_bias_c1 * 0.0
        datagalflag = data["weak_lensing_flag"]

        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_derived_weight, ishape_hsm_regauss_derived_rms_e, ishape_hsm_regauss_derived_bias_m, ishape_hsm_regauss_derived_bias_c1, ishape_hsm_regauss_derived_bias_c2, zeroarr, zeroarr]).T

        hdulist = pyfits.open(sourceargs['fits_list'][itern].split("_pdf_%s" % sourceargs['pofz_type'])[0] + ".fits")
        data = hdulist[1].data
        if sourceargs['pofz_type'] == "frankenz":
            pofzflag = (datagalflag == True)
        elif sourceargs['pofz_type'] == "mizuki":
            pofzflag = (datagalflag == True)
        elif sourceargs['pofz_type'] == "mlz":
            pofzflag = (datagalflag == True)
        elif sourceargs['pofz_type'] == "ephor":
            pofzflag = (datagalflag == True)

        sourceargs['fitspofz'] = fitsio.FITS(sourceargs['fits_list'][itern])

        pofz = sourceargs['fitspofz'][1]['P(z)'][:]
        pofz = pofz.reshape((pofz.shape[0], -1,))

        hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/%s_tracts/%s_pz.fits" % (tractreg, tractnum))
        data = hdulist[1].data
        pofzflag = (data["%s_photoz_risk_best" % (sourceargs['pofz_type']) ]<0.5)

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]
        status = 0

    return datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag

    # S15B v2.1 kept here for notes
    if sourceargs['type'] == "hsc-wide" and sourceargs['filetype'] == "fits" and  sourceargs["fullpofz"]:
        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits_list'] = glob.glob("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S15B_v2.1/*_tracts/*_%s.fits" % (sourceargs['pofz_type']))
            
        if itern == len(sourceargs['fits_list']):
            status = 1

        datagal = 0
        status = not (itern<len(sourceargs['fits_list']))
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0., 0.

        hdulist = pyfits.open(sourceargs['fits_list'][itern].split("_pz_pdf_%s" % (sourceargs['pofz_type']))[0] + ".fits")
        data = hdulist[1].data
        imag_cmodel = data["imag_cmodel"]
        iflux_cmodel = data["iflux_cmodel"]
        iflux_cmodel_err = data["iflux_cmodel_err"]
        ishape_hsm_regauss_resolution = data["ishape_hsm_regauss_resolution"]
        ishape_hsm_regauss_e1 = data["ishape_hsm_regauss_e1"]
        ishape_hsm_regauss_e2 = data["ishape_hsm_regauss_e2"]
        ishape_hsm_regauss_sigma = data["ishape_hsm_regauss_sigma"] 
        ira = data["ira"]
        idec = data["idec"]

        datagalflag = (imag_cmodel < 25.0) & (iflux_cmodel/iflux_cmodel_err >=10) & (ishape_hsm_regauss_resolution> 1./3.) & (ishape_hsm_regauss_e1**2 + ishape_hsm_regauss_e2**2 < 4) & (ishape_hsm_regauss_sigma<0.4)
        #datagalflag = (imag_cmodel < 24.5) & (iflux_cmodel/iflux_cmodel_err >=10) & (iflux_cmodel/iflux_cmodel_err <80) & (ishape_hsm_regauss_resolution> 0.3) & (ishape_hsm_regauss_e1**2 + ishape_hsm_regauss_e2**2 < 4) & (ishape_hsm_regauss_sigma<0.4)
        #datagalflag = (imag_cmodel < 27.0)

        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_sigma]).T

        hdulist = pyfits.open(sourceargs['fits_list'][itern].split("_pdf_%s" % sourceargs['pofz_type'])[0] + ".fits")
        data = hdulist[1].data
        if sourceargs['pofz_type'] == "demp":
            z_median = data["%s_photoz_median" % (sourceargs['pofz_type'])]
            z_mode = data["%s_photoz_mode" % (sourceargs['pofz_type'])]
            z_conf = data["%s_photoz_conf_median" % (sourceargs['pofz_type'])]
            pofzflag = (np.absolute(z_mode - z_median)/(1+z_mode) < 0.8) & (z_conf > 0.09)
        elif sourceargs['pofz_type'] == "mizuki":
            red_chi2 = data["%s_reduced_chisq" % (sourceargs['pofz_type'])]
            z_var = data["%s_photoz_variance" % (sourceargs['pofz_type'])]
            pofzflag = (red_chi2<5) & (z_var<0.45)
        elif sourceargs['pofz_type'] == "nnpz":
            photoz_flag = data["%s_photoz_flag" % (sourceargs['pofz_type'])]
            is_clean = data["%s_is_clean" % (sourceargs['pofz_type'])]
            pofzflag = (is_clean==1) & (photoz_flag<=1)
        elif sourceargs['pofz_type'] == "mlz":
            stddev = data["%s_photoz_stddev_mean" % (sourceargs['pofz_type'])]
            conf_mean = data["%s_photoz_conf_mean" % (sourceargs['pofz_type'])]
            pofzflag = (stddev<3) & (conf_mean>0.13)
        elif sourceargs['pofz_type'] == "frrfv0":
            pofzflag = (datagalflag == True)

        sourceargs['fitspofz'] = fitsio.FITS(sourceargs['fits_list'][itern])

        pofz = sourceargs['fitspofz'][1]['P(z)'][:]
        pofz = pofz.reshape((pofz.shape[0], -1,))

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]
        status = 0

    if sourceargs['type'] == "alexie-test" and sourceargs['filetype'] == "ascii" and  not sourceargs["fullpofz"]:

        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['dfchunks'] = pandas.read_csv("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/alexie-test/GAMA15H_sources_test_original2.cat", names=(['ra','dec','z_b','jk','e1','e2','weight']), chunksize=chunksize, delim_whitespace=1)
        try:
            datagal = sourceargs['dfchunks'].next()
        except:
            status = 1
            return 0.0, sourceargs, 0.0, status, 0, 0, 0

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = datagal.ra.size
        status = 0

        pofz = datagal.z_b.values
        datagalflag = (datagal.ra.values==datagal.ra.values)
        pofzflag = (datagal.ra.values==datagal.ra.values)

    if sourceargs['type'] == "hsc-wide" and sourceargs['filetype'] == "fits" and  not sourceargs["fullpofz"]:

        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits_list'] = glob.glob("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S15B_v2.1/*_tracts/*_pz.fits")
            #sourceargs['fits_list'] = glob.glob("/home/surhud//home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S15B_v2.1/WIDE12H_tracts/*_pz.fits")
            
        if itern == len(sourceargs['fits_list']):
            status = 1

        datagal = 0
        status = not (itern<len(sourceargs['fits_list']))
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0, 0

        # First apply the flags
        hdulist = pyfits.open(sourceargs['fits_list'][itern].split("_pz")[0] + ".fits" )
        data = hdulist[1].data
        imag_cmodel = data["imag_cmodel"]
        iflux_cmodel = data["iflux_cmodel"]
        iflux_cmodel_err = data["iflux_cmodel_err"]
        ishape_hsm_regauss_resolution = data["ishape_hsm_regauss_resolution"]
        ishape_hsm_regauss_e1 = data["ishape_hsm_regauss_e1"]
        ishape_hsm_regauss_e2 = data["ishape_hsm_regauss_e2"]
        ishape_hsm_regauss_sigma = data["ishape_hsm_regauss_sigma"] 
        ira = data["ira"]
        idec = data["idec"]

        datagalflag = (imag_cmodel < 25.0) & (iflux_cmodel/iflux_cmodel_err >=10) & (ishape_hsm_regauss_resolution> 0.3) & (ishape_hsm_regauss_e1**2 + ishape_hsm_regauss_e2**2 < 4) & (ishape_hsm_regauss_sigma<0.4)
        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_sigma]).T

        hdulist = pyfits.open(sourceargs['fits_list'][itern])
        data = hdulist[1].data
        z_median = data["%s_photoz_median" % (sourceargs['pofz_type'])]
        if sourceargs['pofz_type'] == "demp":
            z_mode = data["%s_photoz_mode" % (sourceargs['pofz_type'])]
            z_conf = data["%s_photoz_conf_median" % (sourceargs['pofz_type'])]
            pofzflag = (np.absolute(z_mode - z_median)/(1+z_mode) < 0.8) & (z_conf > 0.09)
        elif sourceargs['pofz_type'] == "mizuki":
            red_chi2 = data["%s_reduced_chisq" % (sourceargs['pofz_type'])]
            z_var = data["%s_photoz_variance" % (sourceargs['pofz_type'])]
            pofzflag = (red_chi2<5) & (z_var<0.45)
        elif sourceargs['pofz_type'] == "nnpz":
            photoz_flag = data["%s_photoz_flag" % (sourceargs['pofz_type'])]
            is_clean = data["%s_is_clean" % (sourceargs['pofz_type'])]
            pofzflag = (is_clean==1) & (photoz_flag<=1)
        elif sourceargs['pofz_type'] == "mlz":
            stddev = data["%s_photoz_stddev_mean" % (sourceargs['pofz_type'])]
            conf_mean = data["%s_photoz_conf_mean" % (sourceargs['pofz_type'])]
            pofzflag = (stddev<3) & (conf_mean>0.13)

        pofz = data["%s_photoz_median" % (sourceargs['pofz_type']) ]

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = pofz.size
        status = 0

    if sourceargs['type'] == "hsc-wide-josh" and sourceargs['filetype'] == "fits":
        # Glob up all the tracts
        if sourceargs["fullpofz"]:
            print "P(z) not implemented for Josh's catalog"
            exit(11)
        else:
            pass

        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits'] = fitsio.FITS("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/RF/hsc_s15b_wl_catalog_v0.fits")
            #sourceargs['fits'] = fitsio.FITS("/home/surhud//home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/rf_wide12h_hsc_s15b_wl_catalog_v0.fits")
            sourceargs['nrows'] = sourceargs['fits'][1].read_header()['naxis2']

        datagal = 0
        status = (itern*chunksize>=sourceargs['nrows'])
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0, 0

        wheremask = sourceargs['fits'][1].where("#row>%d && #row<=%d" % (itern*chunksize, (itern+1)*chunksize))
        datagal = sourceargs['fits'][1]['ra','dec','e1','e2','weight'][wheremask]
        if sourceargs["fullpofz"]:
            print "Full p(z) not supported yet"
        else:
            pofz = sourceargs['fits'][1]['z_med'][wheremask]

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]
        forflag = np.arange(Ngal)
        datagalflag = (forflag==forflag)
        pofzflag = (forflag==forflag)
            
    '''

def get_pofz_array(sourceargs):
    hdulist = pyfits.open("/home/amit/Desktop/Covid_backup/github/weaklens_pipeline/DataStore/S16A_v2.0/pz_pdf_bins_%s.fits" % sourceargs['pofz_type'])
    data = hdulist[1].data
    return data["bins"]
