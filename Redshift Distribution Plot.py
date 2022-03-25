import matplotlib.pyplot as plt
import pandas
import numpy as np
import scipy.stats as stats
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D


fraction = 0.01
#loads in final and initial bin labelings 
final = np.load('/home/hamza/Documents/Results/state_desi_subsample.npy')
initial = np.load('/home/hamza/Documents/UG21_ha/Data/desi_subsample/photometric_labels_initial_hamza.npy')
#files needed to slice redshift values for photometric sample
photometric_sample = np.load('/home/hamza/Documents/UG21_ha/Data/desi_subsample/photometric_sample.npy')
accept = np.load('/home/hamza/Documents/UG21_ha/Data/accept_desi.npy')

#loads in redshift data and slices
true_z = np.load('/home/hamza/Documents/UG21_ha/Data/redshift_true.npy')
true_z = true_z[accept]
true_z = true_z[photometric_sample]

colours = ['b','g','r','orange','purple','grey','black']

#plot preamble
plt.figure(figsize=(16, 8), dpi=100)
plt.rcParams.update({'font.size': 15})
plt.xlabel('z')
plt.ylabel('n(z)')
plt.title('Initial vs final distribution of DESI sample')
#x values needed for kde
x=np.linspace(0,3,500)
patches=[]
#for loop used to plot kde of initial and final redshift distributions for each bin
for j in range(7):
    #finds index for current bin in final results
	indices = np.where(final == j)
	#finds true redshift of the bin
    bin_data = true_z[indices]
    #creates kdes and plots them  	
	kde_final = stats.gaussian_kde(bin_data)
	kde_initial = stats.gaussian_kde(true_z[np.where(initial ==j)])
	plt.plot(x,kde_final.evaluate(x),color=colours[j], label = "final bin {}".format(j+1))
	plt.plot(x,kde_initial.evaluate(x),color=colours[j], linestyle = 'dashed', label = "initial bin {}".format(j+1))
	#creates patches for legend
    patches.append(mpatches.Patch(color=colours[j], label='bin {}'.format(j+1)))

patches.append(Line2D([0], [0],linestyle='--', label='initial', color='k'))
patches.append(Line2D([0], [0],linestyle='-', label='final', color='grey'))

plt.vlines((1/3), 0, 3.7, linestyles = 'dashed', color="gray")
plt.vlines((2/3), 0, 3.7, linestyles = 'dashed', color="gray")
plt.vlines(1, 0, 3.7, linestyles = 'dashed', color="gray")
plt.vlines((4/3), 0, 3.7, linestyles = 'dashed', color="gray")
plt.vlines((5/3), 0, 3.7, linestyles = 'dashed', color="gray")
plt.vlines(2, 0, 3.7, linestyles = 'dashed', color="gray")
plt.vlines(3, 0, 3.7, linestyles = 'dashed', color="gray")

plt.legend(handles=patches)
plt.show()
