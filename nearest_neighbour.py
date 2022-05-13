import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from astropy import units as u
import math as m

### Reading in data ###

data = np.loadtxt('Dendrograms/dendro_values.txt') #read in file of flux data
data_5sig = np.loadtxt('Dendrograms/dendro_values_5sig.txt')
areaell = data[:,1]
areaex = data[:,2]
flux = data[:,3]
majsig = data[:,4]
minsig = data[:,5]
posang = data[:,6]
rad = data[:,7]
xcen = data[:,8]
ycen = data[:,9]

xcen_5sig = data_5sig[:,8]
ycen_5sig = data_5sig[:,9]

### Defining parameters of observation ###

pixwidrad = 1.2120342027737623e-07
d = 8400.0

### Making arrays of positions ###

pos = np.array([xcen,ycen])
posT = np.transpose(pos)
pos_5sig = np.array([xcen_5sig,ycen_5sig])
posT_5sig = np.transpose(pos_5sig)

### Nearest neighbour calculation ###

nbrs = NearestNeighbors(n_neighbors=2, algorithm='auto').fit(posT)
distances, indices = nbrs.kneighbors(posT)

nbrs_5sig = NearestNeighbors(n_neighbors=2, algorithm='auto').fit(posT_5sig)
distances_5sig, indices_5sig = nbrs.kneighbors(posT_5sig)

### Calculating and correcting real distances ###

distrad = distances*pixwidrad
distpc = 2*d*np.tan(distrad/2)
distau = distpc*206265
meandist = np.mean(distpc) * u.pc
meandist_au = meandist.to(u.AU)
distau = np.delete(distau,2,axis=0)
distaunew = distau[:,1]
realdist = 2.3*distaunew
proj_corr_dist = realdist*m.sqrt(3/2)
logdist=np.log10(proj_corr_dist)
mean = np.mean(proj_corr_dist)
std = np.std(proj_corr_dist)
logmean = np.mean(logdist)
logstd = np.std(logdist)

### Calculating and correcting real distances (5 sigma sources) ###

distrad_5sig = distances_5sig*pixwidrad
distpc_5sig = 2*d*np.tan(distrad_5sig/2)
distau_5sig = distpc_5sig*206265
meandist_5sig = np.mean(distpc_5sig) * u.pc
meandist_au_5sig = meandist_5sig.to(u.AU)
distau_5sig = np.delete(distau_5sig,2,axis=0)
distaunew_5sig = distau_5sig[:,1]
realdist_5sig = 2.3*distaunew_5sig
proj_corr_dist_5sig = realdist_5sig*m.sqrt(3/2)
logdist_5sig=np.log10(proj_corr_dist_5sig)
mean_5sig = np.mean(proj_corr_dist_5sig)
std_5sig= np.std(proj_corr_dist_5sig)
logmean_5sig = np.mean(logdist_5sig)
logstd_5sig = np.std(logdist_5sig)

### Defining thermal and turbulent Jeans lengths from previous calculation ###

j_low = 4625.2157761
j_up = 12173.66653489
t_low = 36894.7432496
t_up = 68665.54913149

t_low_half = t_low/2
t_up_half = t_up/2
log_j_low=np.log10(j_low)
log_j_up=np.log10(j_up)
log_t_low=np.log10(t_low)
log_t_up=np.log10(t_up)

log_t_low_half = np.log10(t_low_half)
log_t_up_half = np.log10(t_up_half)

### Plotting histogram ###

plt.figure(figsize=(24,16))
plt.ylim(0,17.5)
plt.plot([logmean,logmean], [-1, 18], 'k--', linewidth=4, label='mean')
plt.hist(logdist,bins=20, color='khaki')
plt.hist(logdist,bins=20, histtype='step', color='black', linewidth=2)
plt.hist(logdist_5sig,bins=5, histtype='step', color='black', linewidth=2, hatch='/', alpha=0.75)
plt.axvspan(log_j_low, log_j_up, facecolor='lightgreen', edgecolor='green', alpha=0.3)
plt.axvspan(log_t_low, log_t_up, facecolor='mediumslateblue', edgecolor='blue', alpha=0.3)
plt.annotate('Thermal Jeans \n Fragmentation', (3.875, 14), fontsize=38, horizontalalignment='center')
plt.annotate('Turbulent \n Fragmentation', (4.7, 14), fontsize=38, horizontalalignment='center')
plt.xlabel('log(separation [au])', fontsize=48, labelpad=22)
plt.ylabel('Frequency', fontsize=48, labelpad=22)
plt.xticks(fontsize=44)
plt.yticks(fontsize=44)
plt.show()
plt.savefig('nearestneighbour.png')