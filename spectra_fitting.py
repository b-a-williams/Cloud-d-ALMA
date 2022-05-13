import radio_beam
import pyspeckit
import numpy as np
from astropy.io import fits
import glob
from astropy import units as u
import matplotlib.pyplot as plt

### Original file in Jy/beam ###

data = fits.getdata('H2CO_30_meanspec.fits')
head = fits.getheader('H2CO_30_meanspec.fits')

### Define beam size and rest frequency ###

bmaj            = (head['bmaj']*u.deg).to(u.arcsec)
bmin            = (head['bmin']*u.deg).to(u.arcsec)
v               = 218.2*u.GHz

### Conversion to get brightness temperature (K) ###

beam    = radio_beam.Beam(bmaj, bmin)

### Creating empty arrays ###

sigma = np.zeros(shape=[94])
centralv = np.zeros(shape=[94])
amplitude = np.zeros(shape=[94])
sigmaerr = np.zeros(shape=[94])
amperr = np.zeros(shape=[94])
velerr = np.zeros(shape=[94])
m = 0

### Loop to fit Gaussian distribution to each extracted spectrum ###

for filename in glob.iglob('H2CO_*''*meanspec.fits'): 
    dat = fits.getdata(filename)
    if len(dat[np.isnan(dat)==True]) > 0:
        print(filename + ' has no data.')
    else:
        sp = pyspeckit.Spectrum(filename, xarrkwargs={'unit':'km/s'})
        sp.xarr.convert_to_unit('km/s')
        sp.unit = "K"
        sp.data = (sp.data*(beam.jtok(v).value))
        d = sp.data
        xaxis = sp.xarr.as_unit('km/s').value
        amp = d.max()
        cenv = xaxis[np.where(d==amp)[0][0]]
        FWHM = ((d[d>0].sum() / amp / (2*np.pi)))**0.5
        guess = [amp, cenv, FWHM]
        sp.specfit(fittype = 'gaussian', Interactive = False, color = 'red', guesses=guess)
        sp.plotter(linestyle='-')
        sp.specfit.plot_fit()
        sp.plotter.savefig("Gaussians/" + filename.strip('.fits') + "_fit.png", overwrite=True)
        ampl = sp.specfit.parinfo['AMPLITUDE0'].value
        amp_err = sp.specfit.parinfo['AMPLITUDE0'].error
        vel = sp.specfit.parinfo['SHIFT0'].value
        vel_err = sp.specfit.parinfo['SHIFT0'].error
        sig = sp.specfit.parinfo['WIDTH0'].value
        sig_err = sp.specfit.parinfo['WIDTH0'].error
        sigma[m] = sig
        centralv[m] = vel
        amplitude[m] = ampl
        sigmaerr[m] = sig_err
        velerr[m] = vel_err
        amperr[m] = amp_err
        m = m+1
        plt.close()
        
### Save values and errors ###
        
np.savetxt('H2CO_sigmas.txt', sigma)
np.savetxt('H2CO_centralvs.txt', centralv)
np.savetxt('H2CO_amplitudes.txt', amplitude)
np.savetxt('H2CO_sigmaerrs.txt', sigmaerr)
np.savetxt('H2CO_centralverrs.txt', velerr)
np.savetxt('H2CO_amplitudeerrs.txt', amperr)