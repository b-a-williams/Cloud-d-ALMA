from __future__ import division
import numpy as np
from astropy import wcs
from astropy import units as u
from astropy.io import fits
from astrodendro import Dendrogram, pp_catalog

### Computing dendrogram as before to subtract background from ###

rms = 1.6e-5 # rms from casaviewer
threesig = 3*rms # 3 sigma sensitivity
file   = 'FITSFiles/d_feather_with_1arcsec.fits'
hdu    = fits.open(file)[0]
dendro = Dendrogram.compute(np.squeeze(hdu.data), min_value=threesig, min_delta=rms, min_npix=77)
dend   = dendro.leaves

h                           = fits.getheader(file)
metadata                    = {}
metadata['data_unit']       = u.Jy / u.beam
metadata['spatial_scale']   = (np.abs(h['CDELT2'])) * u.deg
metadata['beam_major']      = h['bmaj'] * u.deg
metadata['beam_minor']      = h['bmin'] * u.deg
metadata['WCS']             = wcs

cat = pp_catalog(dend, metadata)

### Converting and subtracting background ###

"""
Gives the conversion factor from Jy/beam -> Jy.
Taking the sum of leaf values in Jy/beam and dividing by the
flux (in Jy) from the catalogue to recover the conversion factor.
"""
jy_conversion = np.sum(dend[0].values())/cat['flux'][0]

"""
This loops over each leaf and does the following:
- Gets the minimum value and converts to Jy.
- Multiplies this value by the number of pixels in the leaf.
- Takes the different between the total flux in the catalogue and this.
- Equivalent to doing this on a pixel-by-pixel basis for each leaf.
"""
n       = 0
sub_val  = np.zeros(shape=[96])
sub_flux = np.zeros(shape=[96])
bg_sub = np.zeros(shape=[96])


for leaf in dend:
    sub_val[n]  = leaf.vmin / jy_conversion
    sub_flux[n] = sub_val[n] * leaf.get_npix()
    bg_sub[n]   = cat['flux'][n] - sub_flux[n]
    n = n + 1

np.savetxt('Dendrograms/bs_fluxes.txt', bg_sub)