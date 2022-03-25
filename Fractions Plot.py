import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#load in initial and final states
final = np.load('/home/hamza/Documents/Results/state_desi_subsample.npy')
initial = np.load('/home/hamza/Documents/UG21_ha/Data/desi_subsample/photometric_labels_initial_hamza.npy')
#load in list of accepted galaxies and photometric sample indices
accept = np.load('/home/hamza/Documents/UG21_ha/Data/accept_desi.npy')
phot_sample = np.load('/home/hamza/Documents/UG21_ha/Data/desi_subsample/photometric_sample.npy')

#load in true redshift data
#true_z = pd.read_csv('/home/hamza/Documents/Data/7_bin_data/catalogue.bz2',compression = 'bz2',usecols=['true_z'])
true_z = np.load('/home/hamza/Documents/UG21_ha/Data/redshift_true.npy')
true_z = true_z[accept]
true_z = true_z[phot_sample]
#true_z = np.asarray(true_z['true_z'])
nbins=6
#create bin boundaries
bins  = []
for i in range(nbins+1):
        bins.append(i/(nbins/2))

fraction_in_final = np.empty((nbins+1))
fraction_out_final = np.empty((nbins+1))
fraction_in_initial = np.empty((nbins+1))
fraction_out_initial = np.empty((nbins+1))

for i in range(nbins+1):
    if i == nbins: #checks whether final bin or not
        final_indices = np.where(final == i) #finds indices of current bin in loop for final state
        initial_indices = np.where(initial == i) #finds indices of current bin in loop for initial state
        bin_data_final = true_z[final_indices] #finds redshift data for final state for current bin in loop
        bin_data_initial = true_z[initial_indices] #finds redshift data for initial state for current bin in loop
        fraction_in_final[i] = len(np.where(bin_data_final>=bins[i])[0])/len(final_indices[0]) #calculates fraction of final state inside bin boundary
        fraction_in_initial[i] = len(np.where(bin_data_initial>=bins[i])[0])/len(initial_indices[0]) #calculates fraction of initial state inside bin boundary
        fraction_out_final[i] = len(np.where(bin_data_final<bins[i])[0])/len(final_indices[0]) #calculates fraction of final state outside bin boundary
        fraction_out_initial[i] = len(np.where(bin_data_initial<bins[i])[0])/len(initial_indices[0]) #calculates fraction of initial state outside bin boundary
    else: #does the same as above but for every other bin
        final_indices = np.where(final == i)
        initial_indices = np.where(initial == i)
        bin_data_final = true_z[final_indices]
        bin_data_initial = true_z[initial_indices]
        fraction_in_final[i] = len(np.where((bin_data_final>=bins[i]) & (bin_data_final<bins[i+1]))[0])/len(final_indices[0])
        fraction_in_initial[i] = len(np.where((bin_data_initial>=bins[i]) & (bin_data_initial<bins[i+1]))[0])/len(initial_indices[0])
        fraction_out_final[i] = len(np.where((bin_data_final<bins[i]) | (bin_data_final>=bins[i+1]))[0])/len(final_indices[0])
        fraction_out_initial[i] = len(np.where((bin_data_initial<bins[i]) | (bin_data_initial>=bins[i+1]))[0])/len(initial_indices[0])

#plots bar graphs for the fractions we just calculated 

plt.figure(figsize=(15,10),dpi=80)
plt.rcParams.update({'font.size': 30})
plt.xlabel('Bins',fontsize=38)
plt.ylabel('Fraction',fontsize=38)
plt.title('Fraction of redshift distribution \n within bin boundaries')
width=0.3
plt.bar(np.arange(1,nbins+2), fraction_in_initial, width=width,edgecolor='black', color='lightblue', label='Initial')
#offsets the bars for the final distribution and places the xtick in the middle of both bars
plt.bar(np.arange(1+width,nbins+1.1+width), fraction_in_final, width=width,edgecolor='black', color='blue', label='Final')
plt.xticks(np.arange(1,nbins+2) + width / 2, (np.arange(1,nbins+2)))
plt.legend()
plt.show()

plt.figure(figsize=(15,10),dpi=80)
plt.rcParams.update({'font.size': 30})
plt.xlabel('Bins',fontsize=38)
plt.ylabel('Fraction',fontsize=38)
plt.title('Fraction of redshift distribution \n outside bin boundaries')
width=0.3
plt.bar(np.arange(1,nbins+2), fraction_out_initial, width=width,edgecolor='black', color='lightblue', label='Initial')
plt.bar(np.arange(1+width,nbins+1.1+width), fraction_out_final, width=width,edgecolor='black', color='blue', label='Final')
plt.xticks(np.arange(1,nbins+2) + width / 2, (np.arange(1,nbins+2)))
plt.legend()
plt.show()
