import numpy as np
import pandas as pd
import treecorr
data_folder = '/home/hamza/Documents/UG21_ha/Data'
cat_folder = '/home/hamza/Documents/Data/7_bin_data'
#fraction of reference to photometric galaxies
fraction = 0.01
#loads in accepted galaxies indicies and reference/photometric sample indices
accept = np.load(data_folder + '/accept_desi.npy')
#loads in redshift_true for reference labels, photoz_mode for photometric labels and AllReferences for sample creation
cat = pd.read_csv(cat_folder + '/full_catalogue.bz2', compression='bz2',usecols=['redshift_true','photoz_mode','AllReferences'])
#slices catalogue to accepted galaxies only
cat = cat.iloc[accept]

nbins = 6
#creates bin boundaries
bins  = []
for i in range(nbins+1):
        bins.append(i/(nbins/2))


#creates an array with indices for desi reference and photometric samples
desi_sample=np.where(desi_cat['AllReferences']== True)[0]
phot_sample=np.where(desi_cat['AllReferences']== False)[0]
#saves these arrays
np.save(data_folder + '/reference_sample_%sperc.npy'%fraction,desi_sample)
np.save(data_folder + '/photometric_sample_%sperc.npy'%fraction,phot_sample)

#slices catalogue into photometric and reference samples
phot = cat.iloc[phot_sample]
ref = cat.iloc[desi_sample]
#creates photometric sample labels and reference sample labels
phot_labels=np.digitize(np.asarray(phot['photoz_mode']),bins,right=False)-1
ref_labels=np.digitize(np.asarray(ref['true_z']),bins,right=False)-1
#saves labels
np.save(data_folder + '/desi_subsample/photometric_labels_initial_hamza.npy',phot_labels)
np.save(data_folder + '/reference_labels_hamza.npy',ref_labels)

