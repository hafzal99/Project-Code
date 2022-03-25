import numpy as np
import matplotlib.pyplot as plt


data_folder = '/home/hamza/Documents/UG21_ha/Data/'

nbins = 6
fraction = 0.01
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,5))
for i in range(nbins+1):
	#loads in randoms and plots them
	ref_data = np.load(data_folder + 'randoms/ref_bin%s_%sperc.npy'%(i,fraction))
	phot_data = np.load(data_folder + 'desi_subsample/randoms/phot_bin%s_%sperc.npy'%(i,fraction))
	ax1.scatter(phot_data[0],phot_data[1], s=0.1)
	ax1.set_xlabel('RA (degrees)')
	ax1.set_ylabel('Dec (degrees)')
	ax1.set_title('Photometric Random sample')

	ax2.scatter(ref_data[0],ref_data[1], s=0.1)
	ax2.set_xlabel('RA (degrees)')
	ax2.set_ylabel('Dec (degrees)')
	ax2.set_title('Reference Random Sample')

plt.legend()
plt.show()

