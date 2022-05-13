import pyfits
import numpy as np

### FOR EXECUTION IN CASA ###

d = pyfits.getdata('cloudd-12m7msd-newcasa.fits') #get data from single dish data
h = pyfits.getheader('cloudd-12m7msd-newcasa.fits') #get header from this data

conversion_factor = (225.0/259.0)**3.75 #frequency of our ALMA data/frequency of single dish data

d = d*conversion_factor #scale data

pyfits.writeto('cloudd-12m7msd-newcasa_scaled.fits', d, h) #write scaled data to new FITS image

importfits(fitsimage='cloudd-12m7msd-newcasa_scaled.fits', imagename='cloudd-12m7msd-newcasa_scaled.image')  #import scaled image to CASA

feather(imagename='d_feather_with_1arcsec.image', lowres='cloudd-12m7msd-newcasa_scaled.image', highres='d_TEST.image.fits.galactic') #feather our data with scaled image

exportfits(imagename='d_feather_with_1arcsec.image', fitsimage='d_feather_with_1arcsec.fits') #export new feathered image as a FITS image
