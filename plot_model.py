#!/usr/bin/env python

# Usage: plot_model.py <pattern>
# 
# where <pattern> is the -name setting from a wsclean run.
# A 'modconv' image will be produced which is the image file minus the residual.
# This will then be overlaid as contours over the MFS image to visually vet the 
# model for the purposes of selfcal.


import matplotlib
matplotlib.use('Agg')
import numpy
import aplpy
import pylab
import glob
import sys
from astropy.io import fits
from shutil import copyfile
from optparse import OptionParser


def getImage(fitsfile):
        input_hdu = fits.open(fitsfile)[0]
        if len(input_hdu.data.shape) == 2:
                image = numpy.array(input_hdu.data[:,:])
        elif len(input_hdu.data.shape) == 3:
                image = numpy.array(input_hdu.data[0,:,:])
        else:
                image = numpy.array(input_hdu.data[0,0,:,:])
        return image


def flushFits(newimage,fitsfile):
        f = fits.open(fitsfile,mode='update')
        input_hdu = f[0]
        if len(input_hdu.data.shape) == 2:
                input_hdu.data[:,:] = newimage
        elif len(input_hdu.data.shape) == 3:
                input_hdu.data[0,:,:] = newimage
        else:
                input_hdu.data[0,0,:,:] = newimage
        f.flush()


parser = OptionParser(usage='%prog [options] prefix')
parser.add_option('-l','--level',dest='level',help='Contour level for model in map units (default = 1e-4)',default=1e-4)
parser.add_option('-m','--map',dest='mycmap',help='Colour map for image (default = Greys_r)',default='Greys_r')
parser.add_option('-c','--col',dest='mycol',help='Colour for contours (default = magenta)',default='magenta')
parser.add_option('--pmin',dest='pixmin',help='Minimum for colourscale in map units (default = -1e-4)',default=-1e-4)
parser.add_option('--pmax',dest='pixmax',help='Maximum for colourscale in map units (default = 5e-4)',default=5e-4)
parser.add_option('--pngname',dest='pngname',help='Name of output PNG file (default = something sensible)',default='')
(options,args) = parser.parse_args()


level = float(options.level)
mycmap = options.mycmap
mycol = options.mycol
pixmin = float(options.pixmin)
pixmax = float(options.pixmax)


if len(args) != 1:
    print 'Please specify an image set to plot'
    sys.exit(-1)
else:
    pattern = args[0]

resids = glob.glob('*'+pattern+'*MFS-residual.fits')
first_resids = glob.glob('*'+pattern+'*-first-residual.fits')
fitslist = sorted((f for f in resids if f not in first_resids))

for resid in fitslist:
    img = resid.replace('residual','image')
    modconv = resid.replace('residual','modconv')
    copyfile(img,modconv)
    img_data = getImage(img)
    resid_data = getImage(resid)
    modconv_data = img_data - resid_data
    flushFits(modconv_data,modconv)
    if pngname == '':
            pngname = 'plot_'+img+'.png'
    fig = pylab.figure(figsize=(64,64))
    f1 = aplpy.FITSFigure(img,slices=[0,0],figure=fig,subplot=[0.02,0.02,0.96,0.96])
    f1.show_colorscale(interpolation='none',cmap=mycmap,stretch='linear',vmin=pixmin,vmax=pixmax)
    f1.show_contour(modconv,slices=[0,0],colors=mycol,linewidths=[1.0],alpha=0.7,levels=[level])
    f1.axis_labels.hide()
    f1.tick_labels.hide()
    f1.ticks.hide()
    fig.savefig(pngname,bbox_inches='tight')
    f1.close()


