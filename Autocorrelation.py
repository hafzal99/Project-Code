import numpy as np
import pandas as pd
import treecorr
import gc
data_folder = '/share/splinter/ug21_ha/Data'
cat_folder = '/share/splinter/stolzner/7_bin_data'
#loads in indices of accepted galaxies and the catalogue ra and dec coords, it then slices the catalogue so that we only have accepted galaxies
accept = np.load(cat_folder + '/accept.npy')
cat = pd.read_csv(cat_folder + '/catalogue.bz2', compression='bz2',usecols=['ra','dec'])
cat = cat.iloc[accept]
#loads in reference sample and reference bin labels
ref_sample = np.load(cat_folder + '/reference_sample_0.1perc.npy')
ref_labels = np.load(cat_folder + '/reference_labels.npy')
#sets ra to a variable and slices it to be the ra of the reference galaxies 
ra_ref = np.asarray(cat['ra'])
ra_ref = ra_ref[ref_sample]
#sets dec to a variable and slices it to be the dec of the reference galaxies
dec_ref = np.asarray(cat['dec'])
dec_ref = dec_ref[ref_sample]


fraction = 0.1
nbins = 6
#parameters for correlations
ang_bin = 1
min_sep = 0.01
max_sep = 0.1
sep_units = 'degrees'

#sets arrays for correlations
xi = np.empty((nbins+1))
varxi = np.empty((nbins+1))
rd = np.empty((nbins+1), dtype = treecorr.NNCorrelation)
rr = np.empty((nbins+1), dtype = treecorr.NNCorrelation)
dd = np.empty((nbins+1), dtype = treecorr.NNCorrelation)
dr = np.empty((nbins+1), dtype = treecorr.NNCorrelation)
for i in range(nbins+1):
    #finds ra and dec coords for galaxies in current bin
	ra_data_ref = ra_ref[np.where(ref_labels==i)]
	dec_data_ref = dec_ref[np.where(ref_labels==i)]
	#loads in randoms
    ra_rand_ref, dec_rand_ref = np.load(data_folder + '/randoms_for_auto_20/ref_bin%s_%sperc.npy'%(i,fraction))
	#sets up treecorr catalogues for randoms and data
    rand_ref = treecorr.Catalog(ra = ra_rand_ref, dec = dec_rand_ref, ra_units = 'degrees', dec_units = 'degrees')
	data_ref = treecorr.Catalog(ra = ra_data_ref, dec = dec_data_ref, ra_units = 'degrees', dec_units = 'degrees')
	#calculates the rr correlation for reference bin i with itself
    rr = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=ang_bin,sep_units=sep_units)
	rr.process(rand_ref)
	#calculates the dr correlation for reference bin i with itself
	dr = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=ang_bin,sep_units=sep_units)
	dr.process(data_ref,rand_ref,num_threads = 16)
	#calculates the rd correlation for reference bin i with itself
	rd = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=ang_bin,sep_units=sep_units)
	rd.process(rand_ref,data_ref, num_threads = 16)
	#calculates the dd correlation for reference bin i with itself
	dd = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=ang_bin,sep_units=sep_units)
	dd.process(data_ref,num_threads = 16)
    #evaluates the Landy-Szalay estimator 
	xi[i],varxi[i]=dd.calculateXi(rr, dr, rd)

#saves the autocorrelation and its variance
xi=np.asarray(xi)
varxi = np.asarray(varxi)
np.save(data_folder + '/xi_ref_auto_7_check.npy', xi)
np.save(data_folder + '/varxi_ref_auto_7_check.npy',varxi)
