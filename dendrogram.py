import aplpy
import pyfits
import numpy as np
from astropy import wcs
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.io import fits
from astrodendro import Dendrogram, pp_catalog
from astrodendro.analysis import PPStatistic
from astropy import units as u
from astropy import wcs
from astropy.wcs import WCS

### Specify FITS file and get header ###

file = 'FITSFiles/d_feather_with_1arcsec.fits'
h = pyfits.getheader(file)
data = fits.getdata(file)

mywcs = wcs.WCS(fits.getheader(file))

### Specify beam and pixel size to calculate Npix ###

bmaj = 0.13997140526772
bmin = 0.1094520315528
pixwidrad = 1.2120342027737623e-07
pixwid = 1.2120342027737623e-07*206265
apix = pixwid**2
npix = (2*np.pi*bmaj*bmin)/(8*np.log(2)*apix)

### Specify rms values for dendrogram parameters ###

rms = 1.6e-5
threesig = 3*rms
fivesig = 5*rms

### Compute dendrograms ###

hdu = fits.open(file)[0]
d = Dendrogram.compute(np.squeeze(hdu.data), wcs=mywcs, min_value=threesig, min_delta=rms, min_npix=77)
d.save_to('dendrogram_leaves.fits')

d2 = Dendrogram.compute(np.squeeze(hdu.data), wcs=mywcs, min_value=fivesig, min_delta=rms, min_npix=77)
d2.save_to('dendrogram_leaves_5sig.fits')

### Setting leaves as arrays and creating FITS HDU objects to contain masks ###

mask = np.zeros(hdu.data.shape, dtype=bool)
for leaf in d.leaves:
    mask = mask | leaf.get_mask()

mask_hdu = fits.PrimaryHDU(mask.astype('short'), hdu.header)

mask2 = np.zeros(hdu.data.shape, dtype=bool)
for leaf in d2.leaves:
    mask2 = mask2 | leaf.get_mask()

mask2_hdu = fits.PrimaryHDU(mask2.astype('short'), hdu.header)

### Create a figure of the input file, with dendrogram leaves overlaid as contours. ###

fig = plt.figure(figsize=(24, 20))

f=aplpy.FITSFigure(file, figure=fig)
f.show_colorscale(cmap='Blues')
f.tick_labels.set_xformat('d.ddd')
f.tick_labels.set_yformat('d.ddd')
f.ticks.set_color('black')
f.axis_labels.set_font(size=48)
f.tick_labels.set_font(size=40)
f.axis_labels.set_xpad(8)
f.axis_labels.set_ypad(8)
f.show_contour(mask_hdu, colors='red', linewidths=1)
f.show_contour(mask2_hdu, colors='yellow', linewidths=1)    # Dendrogram leaf contours.
f.set_nan_color('white')

fig.canvas.draw()
fig.savefig('Dendrograms/dendro_figure.png')

### Create a catalogue of the dendrogram leaves and their properties ###

metadata = {}
metadata['data_unit'] = u.Jy / u.beam
metadata['spatial_scale'] = (np.abs(h['CDELT2'])) * u.deg
metadata['beam_major'] = h['bmaj'] * u.deg
metadata['beam_minor'] = h['bmin'] * u.deg
metadata['WCS'] = wcs
cat = pp_catalog(d.leaves, metadata)
ascii.write(cat, 'Dendrograms/dendro_values.txt', format='no_header',overwrite=True)