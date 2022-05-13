#!/usr/bin/env python3

# I haven't tested this script. If it doesn't work, try removing shebang & running in python2.7.
# It's entirely possible I've made a stupid mistake(s) somewhere too.

# Import stuff ...
from __future__ import division
import astrodendro
from spectral_cube import SpectralCube

### Read in dendrogram file from previously ###
dend = astrodendro.Dendrogram.load_from('dendrogram_leaves.fits')
leaves = dend.leaves

### Read in spectral cube - can be swapped out/looped ###
cube_file = 'd_sma1_H2CO.image.pbcor.galactic.fits'
cube = SpectralCube.read(cube_file)

### Exctracting spectra for each dendrogram leaf ###

"""
This will basically loop over all leaves in your dendrogram and extract an averaged spectrum for each.
It takes the mask for each leaf and creates a masked sub-cube.
It then takes an averaged spectrum over each masked cube and writes it out in FITS format.
"""

i = 1
for structure_indx in range(len(leaves)):
    structure = leaves[structure_indx]
    leaf_obj_mask = structure.get_mask()
    leaf_inds = structure.indices()
    view = [slice(leaf_inds[0].min(), leaf_inds[0].max()+1),
            slice(leaf_inds[1].min(), leaf_inds[1].max()+1)]
    submask = leaf_obj_mask[view]
    cubeview = [slice(None),
            slice(leaf_inds[0].min(), leaf_inds[0].max()+1),
            slice(leaf_inds[1].min(), leaf_inds[1].max()+1)]

    cropcube = cube[cubeview].with_mask(submask[None,:,:])
    meanspec = cropcube.mean(axis=(1,2))
    assert meanspec.size == cube.shape[0]
    meanspec.write("H2CO_"+str(i)+"_meanspec.fits", overwrite=True)   # This will write out the averaged spectrum as a FITS file. Read this into pyspeckit to plot/fit/etc. spectrum later.
    # cropcube.write("./meanspec/"+str(i)+"_cube.fits")       # Un-comment this only if you want to write out the masked cubes for each individual leaf.

    i=i+1
