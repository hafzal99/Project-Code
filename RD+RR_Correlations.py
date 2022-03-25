import numpy as np
import pandas as pd
import treecorr
data_folder = '/home/hamza/Documents/UG21_ha/Data'
cat_folder = '/home/hamza/Documents/Data/7_bin_data'

#loads in ra and dec values for galaxies and accepted galaxies list to slice the ra and dec arrays
accept = np.load(cat_folder + '/accept.npy')
cat = pd.read_csv(cat_folder + '/catalogue.bz2', compression='bz2',usecols=['ra','dec'])
cat = cat.iloc[accept]

#loads in reference sample labels and reference sample itself to slice the catalogue
ref_labels = np.load(cat_folder + '/reference_labels.npy')
ref_sample = np.load(cat_folder + '/reference_sample_0.1perc.npy')
ref = cat.iloc[ref_sample]


fraction = 0.1
nbins = 6
#parameters for correlations
ang_bin = 1
min_sep = 0.01
max_sep = 0.1
sep_units = 'degrees'

#initialises arrays for catalogues
ref_cat = np.empty((nbins+1),dtype=treecorr.Catalog)
rand_ref_cat = np.empty((nbins+1),dtype=treecorr.Catalog)
rand_phot_cat = np.empty((nbins+1),dtype=treecorr.Catalog)



for i in range(nbins+1):
    #loads in randoms and creates a catalogue for photometric and reference bins
	rand_phot = np.load(cat_folder + '/new_randoms/phot_bin%s_%sperc.npy'%(i,fraction))
	rand_phot_cat[i] = treecorr.Catalog(ra=rand_phot[0], dec=rand_phot[1], ra_units='degrees', dec_units='degrees')
    
	rand_ref = np.load(cat_folder + '/new_randoms/ref_bin%s_%sperc.npy'%(i,fraction))
	rand_ref_cat[i] = treecorr.Catalog(ra=rand_ref[0], dec=rand_ref[1], ra_units='degrees', dec_units='degrees')
	#finds indices where the galaxies are in this current bin
    ref_bin = ref.iloc[np.where(ref_labels==i)]
    #creates catalogue using current bins galaxies ra and dec values
	ref_cat[i] = treecorr.Catalog(ra=ref_bin['ra'], dec=ref_bin['dec'], ra_units='degrees', dec_units='degrees')
	


for i in range(nbins+1):
		for j in range(nbins+1):
            #calculates rr correlation and writes to disk
			rr = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=ang_bin,sep_units=sep_units)
			rr.process(rand_phot_cat[i],rand_ref_cat[j],num_threads=16)
			rr.write(data_folder + '/rr_correlations_7/rr_phot%s_ref%s_theta%st%s.fits'%(i,j,min_sep,max_sep))
            #calculates rd correlation and writes to disk
			rd = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=nbins,sep_units=sep_units)
			rd.process(rand_phot_cat[i],ref_cat[j],num_threads=16)
			rd.write(data_folder + '/rd_correlations_7/rd_phot%s_ref%s_theta%st%s.fits'%(i,j,min_sep,max_sep))
