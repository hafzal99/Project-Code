import numpy as np
import pandas as pd
import healpy as hp
data_folder = '/share/splinter/ug21_ha/Data'
cat_folder = '/share/splinter/stolzner/7_bin_data'
#loads in accepted galaxies indices and cataloguee
accept = np.load(cat_folder + '/accept.npy')
cat = pd.read_csv(cat_folder + '/full_catalogue.bz2', compression='bz2',usecols=['ra_true','dec_true'])
ref_sample = np.load(data_folder + '/full_cat_ideal_sample/ref_sample.npy')
phot_sample = np.load(data_folder + '/full_cat_ideal_sample/phot_sample.npy')
#loads in reference and photometric bin labels
ref_labels = np.load(data_folder + '/reference_labels.npy')
phot_labels = np.load(data_folder + '/photometric_labels_initial.npy')
#slices catalogue into reference and photometric sample
desi_ref = cat.iloc[ref_sample]
phot = cat.iloc[phot_sample]

fraction = 0.1
nbins = 6 
#healpix pixels of the desi catalogue
pix_ind = [8786, 8787, 8788, 8789, 8790, 8791, 8792, 8793, 8794, 8913, 8914, 8915, 8916, 8917, 8918, 8919, 8920, 8921, 9042, 9043, 9044, 9045, 9046, 9047, 9048, 9049, 9050, 9169, 9170, 9171, 9172, 9173, 9174, 9175, 9176, 9177, 9178, 9298, 9299, 9300, 9301, 9302, 9303, 9304, 9305, 9306, 9425, 9426, 9427, 9428, 9429, 9430, 9431, 9432, 9433, 9434, 9554, 9555, 9556, 9557, 9558, 9559, 9560, 9561, 9562, 9681, 9682, 9683, 9684, 9685, 9686, 9687, 9688, 9689, 9690, 9810, 9811, 9812, 9813, 9814, 9815, 9816, 9817, 9818, 9937, 9938, 9939, 9940, 9941, 9942, 9943, 9944, 9945, 9946, 10066, 10067, 10068, 10069, 10070, 10071, 10072, 10073, 10074, 10193, 10194, 10195, 10196, 10197, 10198, 10199, 10200, 10201, 10202, 10321, 10322, 10323, 10324, 10325, 10326, 10327, 10328, 10329, 10444, 10445, 10446, 10447, 10448, 10449, 10450,10451, 10452]
#a parameter to change the amount of points in the randoms, this is a very bad way to generate these randoms but it seems that if you want m points, then you want to set n=2m
n=20


for i in range(nbins+1):
    #finds galaxies in current bin in the catalogue
    ref_bin = cat.iloc[np.where(ref_labels==i)]
	#maximum and minimum ra and dec	
	ra_min = np.min(ref_bin.ra_true)
	ra_max = np.max(ref_bin.ra_true)
	dec_min = np.min(ref_bin.dec_true)
	dec_max = np.max(ref_bin.dec_true)
    #creates a uniform distribution of ra and dec coordinates between the minimum and maximum ra and dec coords
	rand_ra = np.random.uniform(ra_min, ra_max,int(n*ref_bin.shape[0]))
	rand_sindec = np.random.uniform(np.sin(dec_min*(np.pi/180)),np.sin(dec_max*(np.pi/180)),int(n*ref_bin.shape[0]))
	rand_dec = (180/np.pi)*np.arcsin(rand_sindec)
    #finds the healpix pixels of the ra and dec coords
	rand_pix_ind = hp.ang2pix(32,rand_ra,rand_dec, lonlat = True)
	#checks whether the healpix pixels of the generated ra and dec coords are equal to any of the healpix pixels for the catalogue
	rand_ra = rand_ra[np.where(np.isin(rand_pix_ind,pix_ind) == True)]
	rand_dec = rand_dec[np.where(np.isin(rand_pix_ind,pix_ind)== True)]
	#stacks the ra and dec coords into one array
    rand = np.vstack((rand_ra,rand_dec))
    #saves the randoms to disk
	np.save(data_folder + '/randoms/ref_bin%s_%sperc.npy'%(i,fraction),rand)

#does the same thing as above but for the photometric bins
for i in range(nbins+1):
        
        
        ra_min = np.min(phot_bin.ra_true)
        ra_max = np.max(phot_bin.ra_true)
        dec_min = np.min(phot_bin.dec_true)
        dec_max = np.max(phot_bin.dec_true)
        

        rand_ra = np.random.uniform(ra_min, ra_max,int(n*phot_bin.shape[0]))
        rand_sindec = np.random.uniform(np.sin(dec_min*(np.pi/180)), np.sin(dec_max*(np.pi/180)),int(n*phot_bin.shape[0]))
        rand_dec = (180/np.pi)*np.arcsin(rand_sindec)
        rand_pix_ind = hp.ang2pix(32,rand_ra,rand_dec, lonlat = True)

        rand_ra = rand_ra[np.where(np.isin(rand_pix_ind,pix_ind)==True)]
        rand_dec = rand_dec[np.where(np.isin(rand_pix_ind,pix_ind)==True)]
        rand = np.vstack((rand_ra,rand_dec))

        np.save(data_folder + '/randoms/phot_bin%s_%sperc.npy'%(i,fraction),rand)

