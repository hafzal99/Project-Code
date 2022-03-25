import numpy as np
import pandas as pd
import treecorr
import scipy.stats as stats
import matplotlib.pyplot as plt
import gc
data_folder = '/share/splinter/ug21_ha/Data'
cat_folder = '/share/splinter/stolzner/7_bin_data'
fraction = 0.01
#laods in files to slice redshift and ra/dec values for the samples
accept = np.load(data_folder + '/accept_desi.npy')
ref_sample = np.load(data_folder + '/reference_sample_%sperc.npy'%fraction)
photo_sample = np.load(data_folder + '/desi_subsample/photometric_sample.npy')
#loads in redshift and ra/dec values and slices them 
true_z = np.load(data_folder + '/redshift_true.npy')

true_z = true_z[accept]
true_z = true_z[photo_sample]

ra_true = np.load(data_folder + '/ra_true.npy')
dec_true = np.load(data_folder + '/dec_true.npy')

ra_true = ra_true[accept]
dec_true = dec_true[accept]

ra_ref = ra_true[ref_sample]
dec_ref = dec_true[ref_sample]

ra_phot = ra_true[photo_sample]
dec_phot = dec_true[photo_sample]


#loads in bin labels
phot_labels = np.load(data_folder + '/desi_subsample/photometric_labels_initial_hamza.npy')
ref_labels = np.load(data_folder + '/reference_labels.npy')

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

fig,axes = plt.subplots(nrows=4,ncols=2,figsize=(14.4,24))
x_kde = np.linspace(0,3,1000)
x = np.linspace(1/nbins,2-(1/nbins),nbins)
x = np.append(x,2.5)
#loads in autocorrelation
xi_auto = np.load(data_folder + '/desi_subsample/xi_ref_auto_desi_subsample.npy')
varxi_auto = np.load(data_folder + '/desi_subsample/varxi_ref_auto_desi_subsample.npy')
dp = np.empty((nbins+1,nbins+1))

bins = []
for i in range(nbins+1):
        bins.append(i/(nbins/2))

for i in range(nbins+1):
    #splits photometric sample into bins
    ra_phot_bin = ra_phot[np.where(phot_labels==i)]
    dec_phot_bin = dec_phot[np.where(phot_labels==i)]
    #photometric data catalog
    data_phot[i] = treecorr.Catalog(ra=ra_phot_bin, dec = dec_phot_bin, ra_units = 'degrees', dec_units = 'degrees')

    #splits reference sample into bins
    ra_ref_bin =ra_ref[np.where(ref_labels==i)]
    dec_ref_bin =dec_ref[np.where(ref_labels==i)]
    #loads in reference randoms
    ra_rand_ref, dec_rand_ref = np.load(data_folder + '/randoms/ref_bin%s_%sperc.npy'%(i,fraction))
    #reference randoms catalog
    rand_ref[i] = treecorr.Catalog(ra = ra_rand_ref, dec = dec_rand_ref, ra_units = 'degrees', dec_units = 'degrees')
    #reference data catalog
    data_ref[i] = treecorr.Catalog(ra = ra_ref_bin, dec = dec_ref_bin, ra_units = 'degrees', dec_units = 'degrees')
    
for i in range(nbins+1):
	for j in range(nbins+1):
        #loads in rr correlations for photometric bin i and reference bin j
		rr = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=ang_bin,sep_units=sep_units)
		rr.read(data_folder + '/desi_subsample/rr_correlations/rr_phot%s_ref%s_theta%st%s.fits'%(i,j,min_sep,max_sep))
		#calculates dr correlations for photometric bin i and reference bin j       
		dr = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=ang_bin,sep_units=sep_units)
		dr.process(data_phot[i],rand_ref[j],num_threads=16)
		#loads in rd correlations for photometric bin i and reference bin j
		rd = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=ang_bin,sep_units=sep_units)
		rd.read(data_folder + '/desi_subsample/rd_correlations/rd_phot%s_ref%s_theta%st%s.fits'%(i,j,min_sep,max_sep))
		#calculates dd correlations for photometric bin i and reference bin j
		dd = treecorr.NNCorrelation(min_sep=min_sep,max_sep=max_sep,nbins=ang_bin,sep_units=sep_units)
		dd.process(data_phot[i],data_ref[j],num_threads=16)
		#calculates correlation function using Landy-Szalay estimator
		xi, varxi= dd.calculateXi(rr,dr,rd)
		#calculates cluster z measurements and errors associated with them		
		p[i,j]=xi/(xi_auto[j]**0.5)
		dp[i,j] = p[i,j]*((varxi/xi**2)+(varxi_auto[j]/2*xi_auto[j]**2))**0.5
	if i<nbins:
		if i==0 or i==1:
            #plots first 2 photometric bin cluster z measurements with redshift distribution of the bin itself
			axes[0,i].set_title('{bin1:.2f} < z_phot < {bin2:.2f}'.format(bin1=bins[i],bin2=bins[i+1]),fontsize=15)
			#finds true-z values for galaxies in current bin
            true_z_bin = true_z[np.where(phot_labels==i)]	
			kde = stats.gaussian_kde(true_z_bin)
			kde_y = kde.evaluate(x_kde)
            #calculates factor to multiply cluster-z measurements so that they can be compared with redshift distribution
			factor = np.max(kde_y)/np.max(p[i])
			axes[0,i].errorbar(x,p[i]*factor,yerr=dp[i]*factor,fmt = 'x',color = 'black',label='Cluster-z')
			axes[0,i].plot(x_kde, kde_y,label='n(z) of photometric bin')
			axes[0,i].set_xlabel('z',fontsize=15)
			axes[0,i].set_ylabel('n(z)',fontsize=15)
			axes[0,i].legend()
		if i==2 or i==3:
            #plots photometric bin 3 and 4 cluster z measurements with redshift distribution of the bin itself
			axes[1,i-2].set_title('{bin1:.2f} < z_phot < {bin2:.2f}'.format(bin1=bins[i],bin2=bins[i+1]),fontsize=15)
			true_z_bin = true_z[np.where(phot_labels==i)]
			kde = stats.gaussian_kde(true_z_bin)
			kde_y = kde.evaluate(x_kde)
			factor = np.max(kde_y)/np.max(p[i])
			axes[1,i-2].errorbar(x,p[i]*factor,yerr=dp[i]*factor,fmt = 'x',color = 'black',label='Cluster-z')
			axes[1,i-2].plot(x_kde, kde_y,label='n(z) of photometric bin')
			axes[1,i-2].set_xlabel('z',fontsize=15)
			axes[1,i-2].set_ylabel('n(z)',fontsize=15)
			axes[1,i-2].legend()

		if i==4 or i==5:
            #plots photometric bin 5 and 6 cluster z measurements with redshift distribution of the bin itself
			axes[2,i-4].set_title('{bin1:.2f} < z_phot < {bin2:.2f}'.format(bin1=bins[i],bin2=bins[i+1]),fontsize=15)
			true_z_bin = true_z[np.where(phot_labels==i)]
			kde = stats.gaussian_kde(true_z_bin)
			kde_y = kde.evaluate(x_kde)
			factor = np.max(kde_y)/np.max(p[i])
			axes[2,i-4].errorbar(x,p[i]*factor,yerr=dp[i]*factor,fmt = 'x',color = 'black',label='Cluster-z')
			axes[2,i-4].plot(x_kde, kde_y,label='n(z) of photometric bin')
			axes[2,i-4].set_xlabel('z',fontsize=15)
			axes[2,i-4].set_ylabel('n(z)',fontsize=15)
			axes[2,i-4].legend()
	else:
        #plots photometric bin 7 cluster z measurements with redshift distribution of the bin itself
		axes[3,0].set_title('z_phot > {bin1:.2f}'.format(bin1=bins[i]),fontsize=15)
		true_z_bin = true_z[np.where(phot_labels==i)]
		kde = stats.gaussian_kde(true_z_bin)
		kde_y = kde.evaluate(x_kde)
		factor = np.max(kde_y)/np.max(p[i])
		axes[3,0].errorbar(x,p[i]*factor,yerr=dp[i]*factor,fmt = 'x',color = 'black',label='Cluster-z')
		axes[3,0].plot(x_kde, kde_y,label='n(z) of photometric bin')
		axes[3,0].set_xlabel('z',fontsize=15)
		axes[3,0].set_ylabel('n(z)',fontsize=15)
		axes[3,0].legend()
plt.show()




