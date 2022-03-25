
import numpy as np
import matplotlib.pyplot as plt
#loads in output files from optimisation algorithm for old and new scheme and stacks them 
old_scheme_1 = np.loadtxt('/home/hamza/Documents/Results/og_alg_first10000steps.txt',skiprows=2)
old_scheme_2 = np.loadtxt('/home/hamza/Documents/Results/og_alg_second10000steps.txt',skiprows=2)
old_scheme = np.vstack((old_scheme_1,old_scheme_2))
new_scheme_1 = np.loadtxt('/home/hamza/Documents/Results/new_alg_3700steps.txt',skiprows=2)
new_scheme_2 = np.loadtxt('/home/hamza/Documents/Results/new_alg_second10000steps.txt',skiprows=2)
new_scheme = np.vstack((new_scheme_1,new_scheme_2))
#stacks the time column from both output files
old_scheme_time = np.concatenate((old_scheme_1[:,4],old_scheme_2[:,4]+old_scheme_1[-1,4]))
new_scheme_time = np.concatenate((new_scheme_1[:,4],new_scheme_2[:,4]+new_scheme_1[-1,4]))

#Plots the energy as a function of step
plt.figure(figsize=(12.8,9.6),dpi=80)
plt.rcParams.update({'font.size': 30})
plt.xlabel('Number of Steps')
plt.ylabel('Energy')
plt.title('Energy as a function of steps')

plt.plot(np.arange(old_scheme.shape[0])*100,-old_scheme[:,1],color='blue',label='Old Scheme',linewidth=3)
#line to show when algorithm was restarted for old scheme
plt.vlines(old_scheme_1.shape[0]*100,-np.max(new_scheme[:,1]),-np.min(new_scheme[:,1]),color='blue',linestyle='dashed',label='Old Scheme Restart',linewidth=4)

plt.plot(np.arange(new_scheme.shape[0])*100,-new_scheme[:,1],color='red',label='New Scheme',linewidth=3)
#line to show when algorithm was restarted for new scheme
plt.vlines(new_scheme_1.shape[0]*100,-np.max(new_scheme[:,1]),-np.min(new_scheme[:,1]),color='lightcoral',linestyle='dashed',label='New Scheme Restart',linewidth=4)

plt.legend(prop={'size': 22})
plt.show()

#Plots the energy as a function of time
plt.figure(figsize=(12.8,9.6),dpi=80)
plt.rcParams.update({'font.size': 25})
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('Energy as a function of Time')

plt.plot(old_scheme_time,-old_scheme[:,1],color='blue',label='Old Scheme',linewidth=3)
#line to show when algorithm was restarted for old scheme
plt.vlines(old_scheme_1[-1,4],-np.max(new_scheme[:,1]),-np.min(new_scheme[:,1]),color='blue',linestyle='dashed',label='Old Scheme Restart',linewidth=4)

plt.plot(new_scheme_time,-new_scheme[:,1],color='red',label='New Scheme',linewidth=3)
#line to show when algorithm was restarted for new scheme
plt.vlines(new_scheme_1[-1,4],-np.max(new_scheme[:,1]),-np.min(new_scheme[:,1]),color='lightcoral',linestyle='dashed',label='New Scheme Restart',linewidth=4)

plt.legend(prop={'size': 22})
plt.show()