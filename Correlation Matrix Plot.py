import numpy as np
import pandas as pd
import treecorr
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib
import gc
fraction = 0.01
#files needed to slice catalogue into the samples
accept = np.load('/home/hamza/Documents/UG21_ha/Data/accept_desi.npy')
photometric_sample = np.load('/home/hamza/Documents/UG21_ha/Data/desi_subsample/photometric_sample.npy')
reference_sample = np.load('/home/hamza/Documents/UG21_ha/Data/reference_sample_%sperc.npy'%fraction)
ra = np.load('/home/hamza/Documents/UG21_ha/Data/ra_true.npy')
dec = np.load('/home/hamza/Documents/UG21_ha/Data/dec_true.npy')
#splits sample into photometric and reference sample
ra = ra[accept]
dec = dec[accept]
ra_phot = ra[photometric_sample]
ra_ref = ra[reference_sample]
dec_phot = dec[photometric_sample]
dec_ref = dec[reference_sample]
#loads in results and labels
results = np.load('/home/hamza/Documents/Results/state_desi_subsample.npy')
ref_labels = np.load('/home/hamza/Documents/UG21_ha/Data/reference_labels.npy')
nbins = 6

ang_bin = 1 
min_sep = 0.03
max_sep = 0.3
sep_units = 'degrees'

#correlation matrix
p = np.empty((nbins+1,nbins+1))
#arrays for treecorr catalogs
rand_phot = np.empty((nbins+1), dtype = treecorr.Catalog)
data_phot = np.empty((nbins+1), dtype = treecorr.Catalog)
rand_ref = np.empty((nbins+1), dtype = treecorr.Catalog)
data_ref = np.empty((nbins+1), dtype = treecorr.Catalog)

xi_auto = np.empty((nbins+1))
for i in range(nbins+1):
    #splits photometric sample into bins
    ra_data_phot = ra_phot[np.where(results==i)]
    dec_data_phot = dec_phot[np.where(results==i)]
    #photometric data catalog
    data_phot[i] = treecorr.Catalog(ra = ra_data_phot, dec = dec_data_phot, ra_units = 'degrees', dec_units = 'degrees')

    #splits reference sample into bins
    ra_data_ref = ra_ref[np.where(ref_labels==i)]
    dec_data_ref = dec_ref[np.where(ref_labels==i)]
    
    #loads in reference randoms
    ra_rand_ref, dec_rand_ref = np.load('/home/hamza/Documents/UG21_ha/Data/randoms/ref_bin%s_%sperc.npy'%(i,fraction))
    #reference randoms catalog
    rand_ref[i] = treecorr.Catalog(ra = ra_rand_ref, dec = dec_rand_ref, ra_units = 'degrees', dec_units = 'degrees')
    #reference data catalog
    data_ref[i] = treecorr.Catalog(ra = ra_data_ref, dec = dec_data_ref, ra_units = 'degrees', dec_units = 'degrees')


xi_auto = np.load('/home/hamza/Documents/UG21_ha/Data/desi_subsample/xi_ref_auto_desi_subsample.npy')

for i in range(nbins+1):
	for j in range(nbins+1):
        #laods in rr correlations for photometric bin i with reference bin j
		rr = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=ang_bin,sep_units=sep_units)
		rr.read('/home/hamza/Documents/UG21_ha/Data/desi_subsample/rr_correlations/rr_phot%s_ref%s_theta%st%s.fits'%(i,j,min_sep,max_sep))
        #calculates dr correlations for photometric bin i with reference bin j
		dr = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=ang_bin,sep_units=sep_units)
		dr.process(data_phot[i],rand_ref[j],num_threads=16)
        #loads in rd correlations for photometric bin i with reference bin j
		rd = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=ang_bin,sep_units=sep_units)
		rd.read('/home/hamza/Documents/UG21_ha/Data/desi_subsample/rd_correlations/rd_phot%s_ref%s_theta%st%s.fits'%(i,j,min_sep,max_sep))
        #calculates dd correlations for photometric bin i with reference bin j
		dd = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=ang_bin,sep_units=sep_units)
		dd.process(data_phot[i],data_ref[j],num_threads=16)
        #calculates correlation function using Landy-Szalay estimator
		xi, varxi= dd.calculateXi(rr,dr=dr,rd=rd)
		
		p[i,j]=xi/xi_auto[j]

#plots matrix
fig, ax = plt.subplots(figsize = (8,6),dpi=100)
matplotlib.rcParams.update({'font.size': 15})
mat = ax.matshow(p, cmap='Blues',vmin=np.min(p),vmax=1)
plt.colorbar(mat,ax=ax)
ax.set_xticklabels(['']+list(map(str,[1,2,3,4,5,6,7])),fontsize=18)
ax.set_yticklabels(['']+list(map(str,[1,2,3,4,5,6,7])),fontsize=18)
ax.set_ylabel('Photometric Bin',fontsize=18)
ax.set_xlabel('Reference Bin',fontsize=18)
ax.xaxis.set_ticks_position("bottom")
plt.title('Final correlation matrix for DESI subsample')
print(xi_auto)
for i in range(nbins+1):
    for j in range(nbins+1):
        c = p[j,i]
        if p[j,i]>0.6:
            ax.text(i, j, '{number:.2f}'.format(number=c), va='center', ha='center', color = 'white')
        else:
            ax.text(i, j, '{number:.2f}'.format(number=c), va='center', ha='center', color = 'black')
plt.show()
