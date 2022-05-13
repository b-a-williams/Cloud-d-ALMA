import numpy as np
from matplotlib import pyplot as plt
import math as m
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

### Reading in data from dendrogram calculations ###

data = np.loadtxt('Dendrograms/dendro_values.txt') #read in file of flux data
bs_flux = np.loadtxt('Dendrograms/bs_fluxes.txt')
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

areaex_5sig = data_5sig[:,2]
flux_5sig = data_5sig[:,3]
xcen_5sig = data_5sig[:,8]
ycen_5sig = data_5sig[:,9]

### Regridding pixel positions to Galactic coordinates and printing to check ###

f = fits.open('FITSFiles/d_sma1_ALL_LINEFREE_CONT.image.fits')
w = WCS(f[0].header)
sky = w.pixel_to_world(xcen, ycen)
skyT = np.transpose(sky)
print(skyT)

### Mass calculation ###

T = 20.0 #gas temp in K
d = 8100.0 #distance to galactic centre in pc
beta = 1.75
v = 2.25e+11 #observing frequency in Hz
vghz = v/10e8
k = 0.04*(((v/10e8)/505)**beta)
c = 299792458.0 #speed of light in m/s
L1 = c/v #wavelength in m
L2 = L1*1000.0 #wavelength in mm
msol = 1.989e+30 #solar mass in kg
a = 0.12 #constant from eqn
b = (np.exp(1.439*(1.0/L2)*(10.0/T))) - 1 #first term of eqn
e = a*b*(0.01/k)*((d/100.0)**2)*(L2**3) #conversion factor in msol
mconv = flux*e * u.solMass
m_5sig = flux_5sig*e * u.solMass

bs_mass = bs_flux*e

### Mass sensitivity ###

rms = 1.6e-5 # rms from casaviewer
threesig = 3*rms
threesigminmass = threesig*e*77

### Effective radius ###

radius = np.sqrt(areaex/np.pi)
pcconv = 3.086E+16 #conversion factor for pc to m
radrad = radius*(m.pi/180) #converting radius to radians
radpc = 2*d*np.tan(radrad/2) #converting radius to pc
radau = radpc*206265

radius_5sig = np.sqrt(areaex_5sig/np.pi)
radrad_5sig = radius_5sig*(m.pi/180) #converting radius to radians
radpc_5sig = 2*d*np.tan(radrad_5sig/2) #converting radius to pc
radau_5sig = radpc_5sig*206265

### Density ###

pixwid = 8.333333333333E-06
pixwidrad = pixwid*(m.pi/180) #converting pixel width from degrees to radians
pixwidcm = 2*d*pcconv*np.tan(pixwidrad/2) #converting pixel width to cm
pixareacm = pixwidcm**2

mh2 = 3.35E-27 * u.kg #mass of H2 molecule in kg
mh2msol = mh2.to(u.solMass)
nh2 = mconv/mh2msol
bs_nh2 = bs_mass/mh2msol
nh2_5sig = m_5sig/mh2msol
pcconv = 3.086E+18 #conversion factor of pc to cm
radcm = radpc*pcconv * u.cm
radcm_5sig = radpc_5sig*pcconv * u.cm
dens = nh2/((4/3)*m.pi*(radcm**3)) #calculating densities
bs_dens = bs_nh2/((4/3)*m.pi*(radcm**3))
dens_5sig = nh2_5sig/((4/3)*m.pi*(radcm_5sig**3))
leafno = np.arange(1,97,1)
radaucor = radau/1e3
denscor = dens/1e6
bs_denscor = bs_dens/1e6

### Rounding data down to two decimal places for tables ###

mconv2dp = [round(mass, 2) for mass in mconv.value]
radaucor2dp = [round(rad, 2) for rad in radaucor]
denscor2dp = [round(dens, 2) for dens in denscor.value]
bs_mass2dp = [round(mass, 2) for mass in bs_mass]
bs_denscor2dp = [round(dens, 2) for dens in bs_denscor.value]

### Saving data ###
np.savetxt('Dendrograms/masses.txt', mconv.value)
np.savetxt('Dendrograms/radiipc.txt', radpc)
np.savetxt('Dendrograms/radiiau.txt', radau)
np.savetxt('Dendrograms/densities.txt', dens.value)
np.savetxt('Dendrograms/bs_masses.txt', bs_mass)
np.savetxt('Dendrograms/bs_densities.txt', bs_dens.value)
np.savetxt('Dendrograms/masses_5sig.txt', m_5sig.value)
np.savetxt('Dendrograms/radiipc_5sig.txt', radpc_5sig)
np.savetxt('Dendrograms/radiiau_5sig.txt', radau_5sig)
np.savetxt('Dendrograms/densities_5sig.txt', dens_5sig.value)
np.savetxt("leaf_props.csv",(leafno,skyT,mconv2dp,radaucor2dp,denscor2dp,bs_mass2dp,bs_denscor2dp),delimiter=",",fmt="%s")

### Calculating Salpeter line ###

k1 = 96/(((np.max(mconv.value)**(-1.35))/(-1.35))-((np.min(mconv.value)**(-1.35))/(-1.35)))
Msal = np.linspace(-50, 50, int(1e5))
logMsal = np.log10(Msal)
Nsal = (-1.35*logMsal) + 0.75

### Logging values for plotting ###

logmconv = np.log10(mconv.value)
logbsm = np.log10(bs_mass)
logm5sig = np.log10(m_5sig.value)
mj_low = 0.19859872
mj_up = 1.04543214
log_mj_low = np.log10(mj_low)
log_mj_up = np.log10(mj_up)

### Defining histogram parameters ###

n, bins = np.histogram(logmconv, bins=50)
bin_centres = (bins[:-1] + bins[1:])*0.5
err = np.sqrt(n)
n_bs, bins_bs = np.histogram(logbsm, 50)

### Plotting mass distribution with Salpeter line ###

plt.figure(figsize=(24,16))
plt.xlim(-0.71,0.6)
plt.ylim(0, 17)
plt.plot(logMsal, Nsal, color='blue', linewidth=4)
plt.legend(['Salpeter (1955)'], fontsize=44)
plt.hist(logmconv, bins=20, color='powderblue')
plt.hist(logmconv,bins=20, histtype='step', color='black', linewidth=2)
plt.hist(logm5sig,bins=5, histtype='step', color='black', linewidth=2, hatch='/', alpha=0.75)
plt.axvspan(log_mj_low, log_mj_up, facecolor='pink', edgecolor='red', alpha=0.3)
plt.xlabel('log(M/M${}_\odot$)', fontsize = 48, labelpad=22)
plt.ylabel('dN/d(logM)', fontsize = 48, labelpad=22)
plt.xticks(fontsize=44)
plt.yticks(fontsize=44)
plt.show()
plt.draw()
plt.savefig('massdist.png')