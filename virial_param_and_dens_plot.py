import numpy as np
from astropy import units as u
from astropy import constants as const
import math as m
from matplotlib import pyplot as plt

### Read in data ###

disc = np.loadtxt('disc_points.txt')
discx = disc[:,0]
discy = disc[:,1]

massmsol = np.loadtxt('Dendrograms/masses.txt') * u.solMass
radiipc = np.loadtxt('Dendrograms/radiipc.txt') * u.pc
cenvh2co = np.loadtxt('Spec Fits Files/H2CO_centralvs.txt')
amph2co = np.loadtxt('Spec Fits Files/H2CO_amplitudes.txt')
sigh2cokm = np.loadtxt('Spec Fits Files/H2CO_sigmas.txt') #may not be needed - check
cenverrh2co = np.loadtxt('Spec Fits Files/H2CO_centralverrs.txt')
amperrh2co = np.loadtxt('Spec Fits Files/H2CO_amplitudeerrs.txt')
sigerrh2cokm = np.loadtxt('Spec Fits Files/H2CO_sigmaerrs.txt') #may not be needed - check

### Defining constants ###

G = const.G.value
kB       = const.k_B.value
pc      = const.pc.value
Gamma   = 0.73

### Zero pressure term (standard Virial eq.) ###

n   = 0
Yp  = np.zeros(shape=[5000])
Sig = np.zeros(shape=[5000])
Ypc  = np.zeros(shape=[5000])
Sigc  = np.zeros(shape=[5000])
Pe  = 0

for sigma in np.logspace(-2.5,1.5,num=5000):
    Yp[n]   = (1/3)*((np.pi*Gamma*G*sigma*10) + ((4*(Pe*(kB*1e6)))/(sigma*10)))
    Sig[n]  = sigma
    n       = n+1

### Pressure term included ###

n   = 0
Pe  = [1e5, 1e6, 1e7, 1e8, 1e9]

for i in range(0,len(Pe),1):
    for sigma in np.logspace(-2.5,1.5,num=5000):
        Ypc[n]   = (1/3)*((np.pi*Gamma*G*sigma*10) + ((4*(Pe[i]*(kB*1e6)))/(sigma*10)))
        Sigc[n]  = sigma
        n       = n+1

    n=0
    plt.plot(np.log10(Sigc), np.log10((Ypc*pc)/1e6), color='k', alpha=0.4)

### Removing sources for which we could not extract any spectra ###

indices = (1,2)

masskg = massmsol.to(u.kg)
masskgex = np.delete(masskg, indices, axis=0)
radiim = radiipc.to(u.m)
radiimex = np.delete(radiim, indices, axis=0)

### Calculating velocity dispersions, errors and upper & lower values ###

sigh2co = sigh2cokm*1000
sigerrh2co = sigerrh2cokm*1000

sigh2co_up = sigh2co + sigerrh2co
sigh2co_low = sigh2co - sigerrh2co
sigh2co_up_km = sigh2cokm + sigerrh2cokm
sigh2co_low_km = sigh2cokm - sigerrh2cokm

### Calculating virial parameters and upper & lower values

alphah2co = (((5*radiimex*(sigh2co**2))/(G*masskgex))).value
alphah2co_up = (((5*radiimex*(sigh2co_up**2))/(G*masskgex))).value
alphah2co_low = (((5*radiimex*(sigh2co_low**2))/(G*masskgex))).value     

### Reading in data from dendrogram leaves ###         

data = np.loadtxt('Dendrograms/dendro_values.txt') #read in file of flux data
areaell = data[:,1]
areaex = data[:,2]
flux = data[:,3]
majsig = data[:,4]
minsig = data[:,5]
posang = data[:,6]
rad = data[:,7]
xcen = data[:,8]
ycen = data[:,9]

### Calculating mass of each dendrogram leaf ###

T = 20.0 #gas temp in K
k = 0.01 #opacity of gas in cm2/g
d = 8400.0 #distance to galactic centre in pc
v = 2.3e+11 #observing frequency in Hz
c = 299792458.0 #speed of light in m/s
L1 = c/v #wavelength in m
L2 = L1*1000.0 #wavelength in mm
msol = 1.989e+30 #solar mass in kg
a = 0.12 #constant from eqn
b = (np.exp(1.439*(1.0/L2)*(10.0/T))) - 1 #first term of eqn
e = a*b*((d/100.0)**2)*(L2**3) #conversion factor in msol
mconv = flux*e * u.solMass
mg = mconv.to(u.g)
mgex = np.delete(mg, indices, axis=0)
mconvex = np.delete(mconv, indices, axis=0)

### Calculating aspect ratios of each leaf ###

majsigex = np.delete(majsig, indices, axis=0)
minsigex = np.delete(minsig, indices, axis=0)
asprat = majsigex/minsigex

### Calculating areas of each leaf to calculate surface density ###

radius = np.sqrt(areaex/np.pi)
pcconv = 3.086E+16 #conversion factor for pc to m
radrad = radius*(m.pi/180) #converting radius to radians
radpc = 2*d*np.tan(radrad/2) #converting radius to pc
radpcex = np.delete(radpc, indices, axis=0)
radau = radpc*206265
pcconv = 3.086E+18
radcm = radpc*pcconv
areacm = np.pi*(radcm**2)
areacmex = np.delete(areacm, indices, axis=0)
areapcex = np.pi*(radpcex**2)

Sigma = mgex/areacmex
Sigmamsol = mconvex/areapcex

### Calculating y axis values ###

yaxh2co = sigh2cokm**2/radpcex
yaxh2co_up = (sigh2co_up_km**2)/radpcex
yaxh2co_low = (sigh2co_low_km**2)/radpcex

### Excluding values again ###

ind2 = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,20,21,22,23,24,25,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,48,49,50,51,54,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,80,81,83,84,85,86,87,88,89,90,91,92,93,94)

merr = 2 * mconvex
merrex = np.delete(merr, ind2, axis=0)
areacmexex = np.delete(areacmex, ind2, axis=0)
merrg = merrex.to(u.g)
errSig = merrg/areacmexex
logerrSig = np.log10(errSig.value)

yaxh2co_ex = np.delete(yaxh2co, ind2, axis=0)
yaxh2co_up_ex = np.delete(yaxh2co_up, ind2, axis=0)
yaxh2co_low_ex = np.delete(yaxh2co_low, ind2, axis=0)

### Calculating upper and lower values of y axis ###

yaxh2co_up_err = yaxh2co_up_ex - yaxh2co_ex
yaxh2co_low_err = yaxh2co_ex - yaxh2co_low_ex
log_y_err_up = np.log10(yaxh2co_up_err)
log_y_err_low = np.log10(yaxh2co_low_err)

### Reading in SMA values ###

X = [-0.62087078, -0.47954163, -0.3067933 ,  0.11148861, -0.02655597, 0.08395594]
Y = [1.32944709, 1.32944709, 1.17524341, 1.25647033, 1.15157407, 1.38995607]

### Creating plot ###

n   = 0
Pe  = [1e5, 1e6, 1e7, 1e8, 1e9]

cloud_d_x = np.log10(Sigma[16].value),np.log10(Sigma[19].value),np.log10(Sigma[26].value),np.log10(Sigma[27].value),np.log10(Sigma[47].value),np.log10(Sigma[52].value),np.log10(Sigma[53].value),np.log10(Sigma[55].value),np.log10(Sigma[76].value),np.log10(Sigma[77].value),np.log10(Sigma[78].value),np.log10(Sigma[79].value),np.log10(Sigma[82].value)
cloud_d_y = np.log10(yaxh2co[16]),np.log10(yaxh2co[19]),np.log10(yaxh2co[26]),np.log10(yaxh2co[27]),np.log10(yaxh2co[47]),np.log10(yaxh2co[52]),np.log10(yaxh2co[53]),np.log10(yaxh2co[55]),np.log10(yaxh2co[76]),np.log10(yaxh2co[77]),np.log10(yaxh2co[78]),np.log10(yaxh2co[79]),np.log10(yaxh2co[82])

plt.figure(figsize=(24,16))
plt.errorbar(cloud_d_x, cloud_d_y, yerr=(log_y_err_low,log_y_err_up), xerr=logerrSig, marker='D', color='c', markersize=10, linewidth=0, ecolor='black', elinewidth=1, capsize=6, capthick=1)
plt.plot([-0.62087078,0.11148861,0.08395594], [1.32944709,1.25647033,1.38995607], marker='*', color='g', markersize=10, linewidth=0)
plt.plot(-0.47954163,1.32944709, marker='o', color='b', markersize=10, linewidth=0)
plt.plot(-0.3067933,1.17524341, marker='p', color='r', markersize=10, linewidth=0)
plt.plot(-0.02655597,1.15157407, marker='v', color='m', markersize=10, linewidth=0)
plt.plot(discx, discy, marker='+', color='k', markersize=10, linewidth=0)
plt.plot(np.log10(Sig), np.log10((Yp*pc)/1e6), color='k', ls='--', alpha=0.4, linewidth=2)
for i in range(0,len(Pe),1):
    for sigma in np.logspace(-2.5,1.5,num=5000):
        Ypc[n]   = (1/3)*((np.pi*Gamma*G*sigma*10) + ((4*(Pe[i]*(kB*1e6)))/(sigma*10)))
        Sigc[n]  = sigma
        n       = n+1

    n=0
    plt.plot(np.log10(Sigc), np.log10((Ypc*pc)/1e6), color='k', alpha=0.4, linewidth=2.5)
plt.xlim(-2.5,0.5)
plt.xticks(fontsize=40)
plt.yticks(fontsize=40)
plt.legend(('"d" (ALMA sources)', '"a", "b" & "f" (clouds)*', '"c" (cloud)*', '"d" (cloud)*', '"e" (cloud)*', 'FBK (2011)'), loc='lower right', fontsize=28, ncol=2)
plt.annotate('P$_e/k$ = 10$^5$ K cm$^{-3}$', (-2.45, 0.25), fontsize=30)
plt.annotate('P$_e/k$ = 10$^6$ K cm$^{-3}$', (-2.4, 1.25), fontsize=30)
plt.annotate('P$_e/k$ = 10$^7$ K cm$^{-3}$', (-2, 1.8), fontsize=30)
plt.annotate('P$_e/k$ = 10$^8$ K cm$^{-3}$', (-1.6, 2.4), fontsize=30)
plt.annotate('P$_e/k$ = 10$^9$ K cm$^{-3}$', (-0.9, 2.7), fontsize=30)
plt.xlabel('log $\Sigma$ (g cm$^2$)', fontsize=48)
plt.ylabel('log $\sigma^2$ / R (km$^2$ s$^{-2}$ pc$^{-1}$)', fontsize=48)
plt.savefig('h2co_surface_density.png')

### Excluding values from virial parameters ###

alphah2coex = np.delete(alphah2co, ind2, axis=0)
alphah2co_low_ex = np.delete(alphah2co_low, ind2, axis=0)
alphah2co_up_ex = np.delete(alphah2co_up, ind2, axis=0)

### Calculating errors of virial parameters ###

alphah2co_up_err = alphah2co_up_ex - alphah2coex
alphah2co_low_err = alphah2coex - alphah2co_low_ex
logalphah2co_up_err = np.log10(alphah2co_up_err)
logalphah2co_low_err = np.log10(alphah2co_low_err)

### Plotting aspect ratios vs virial parameters ###

aspratex = np.delete(asprat, ind2, axis=0)
logasprat = np.log10(aspratex)
logalphah2coex = np.log10(alphah2coex)
logalphah2co_low_ex = np.log10(alphah2co_low_ex)
logalphah2co_up_ex = np.log10(alphah2co_up_ex)

plt.figure(figsize=(24,16))
plt.errorbar(logasprat, logalphah2coex, xerr=None, yerr=(logalphah2co_low_err,logalphah2co_up_err), marker='x', markersize=20, linewidth=0, color='black', ecolor='black', elinewidth=1, capsize=6, capthick=1 )
plt.xlabel('log(Major axis / Minor axis)', fontsize = 48, labelpad = 22)
plt.ylabel(r'log($\alpha$)', fontsize = 48, labelpad = 22)
plt.xticks(fontsize=44)
plt.yticks(fontsize=44)
plt.show()

### Writing values to csv file ###

leafno = np.arange(1,95,1)
np.savetxt("virial.csv",(leafno,alphah2co,alphah2co_low,alphah2co_up,sigh2cokm,sigerrh2cokm),delimiter=",",fmt="%f")