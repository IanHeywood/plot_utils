#!/usr/bin/env python

import numpy
import pylab
import matplotlib.colors as colors
import matplotlib.cm as cmx
import pyrap.tables
import glob
import os
import sys
from optparse import OptionParser
from scipy.interpolate import UnivariateSpline


def gi(message):
        print '\033[92m'+message+'\033[0m'


def ri(message):
        print '\033[91m'+message+'\033[0m'

	
# Command line options

parser = OptionParser(usage='%prog [options] gaintable')
parser.add_option('-p','--plot',dest='doplot',help='Set to amp, phase, real or imag (default = amp)',default='amp')
parser.add_option('-c','--corr',dest='corr',help='Correlation product to plot (default = 0)',default=0)
parser.add_option('--ms',dest='msname',help='Corresponding Measurement Set from which to read antenna info (not essential)',default='')
parser.add_option('--showflags',dest='showflags',action='store_false',help='Show flagged points in regular colour, not magenta (default = True)',default=True)
parser.add_option('--y0',dest='y0',help='Minimum value of y-axis',default='')
parser.add_option('--y1',dest='y1',help='Maximum value of y-axis',default='')
parser.add_option('--pngname',dest='pngname',help='Name of output PNG',default='')
(options,args) = parser.parse_args()

doplot = options.doplot
corr = options.corr
msname = options.msname
showflags = options.showflags
y0 = options.y0
y1 = options.y1
pngname = options.pngname


# Options Parsing


if len(args) != 1:
	ri('Please specify a gain table to plot.')
	sys.exit(-1)
else:
	caltab = args[0].rstrip('/')

if doplot not in ['amp','phase','real','imag']:
	ri('Requested plot not valid, must be one of amp, phase, real or imag.')
	sys.exit(-1)

if msname == '':
	automs = glob.glob('*'+caltab.split('.')[0].split('cal_')[-1]+'*.ms')
	if len(automs) > 0 and os.path.isdir(automs[0]):
		gi('Auto-detected '+automs[0])
		msname = automs[0]
		gotms = True
	else:
		ri('Could not autodetect associated MS, antenna names will not be available')
		gotms = False
else:
	gotms = True

if gotms:
	anttab = pyrap.tables.table(msname+'/ANTENNA')
	antnames = anttab.getcol('NAME')
	anttab.done()

# Set y-ranges to some defaults if not specified

if y0 =='':
	if doplot in ['real','amp']:
		ymin = 0.8 
	elif doplot == 'imag':
		ymin = -0.8
	elif doplot == 'phase':
		ymin = -numpy.pi
else:
	ymin = float(y0)

if y1 =='':
	if doplot in ['real','amp']:
		ymax = 1.2 
	elif doplot == 'imag':
		ymax = 0.3
	elif doplot == 'phase':
		ymax = numpy.pi
else:
	ymax = float(y1)

# Set the pngname to a default if not specified

if pngname == '':
	pngname = caltab+'_'+doplot+'_corr'+str(corr)+'.png'

# Open the caltable and gather some basic info

t = pyrap.tables.table(caltab)
spws = numpy.unique(t.getcol('SPECTRAL_WINDOW_ID')).tolist()
ants = numpy.unique(t.getcol('ANTENNA1')).tolist()
scans = numpy.unique(t.getcol('SCAN_NUMBER')).tolist()

# Set colourmap based on number of SPWs

cNorm = colors.Normalize(vmin=0,vmax=len(spws)-1)
mymap = cm = pylab.get_cmap('jet_r')
scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=mymap)

# Open the figure and set the panel grid based on the number of antennas and SPWs

fig = pylab.figure(figsize=(50,50))
c = 1
ncols = len(spws)
#nrows = int(numpy.ceil(float(len(ants))/float(ncols)))
nrows = len(ants)
for ant in ants:
	gi('Antenna '+str(ant))
	for spw in spws:
		ax = fig.add_subplot(nrows,ncols,c,axisbg='#EEEEEE')
		subtab = t.query(query='ANTENNA1=='+str(ant)+' && SPECTRAL_WINDOW_ID=='+str(spw))
		subtimes = subtab.getcol('TIME').tolist()
		subgains = subtab.getcol('CPARAM').tolist()
		subflags = subtab.getcol('FLAG')
		unflaggedgains = numpy.ma.masked_where(subflags==True,subgains)
		flaggedgains = numpy.ma.masked_where(subflags==False,subgains)
		if doplot == 'real':
			unflaggedplot = unflaggedgains.real
			flaggedplot = flaggedgains.real
		elif doplot == 'imag':
			unflaggedplot = unflaggedgains.imag
			flaggedplot = flaggedgains.imag
		elif doplot == 'amp':
			unflaggedplot = numpy.ma.absolute(unflaggedgains)
			flaggedplot = numpy.ma.absolute(flaggedgains)
		elif doplot == 'phase':
			unflaggedplot = numpy.ma.angle(unflaggedgains)
			flaggedplot = numpy.ma.angle(flaggedgains)
		if len(subgains) > 0:
			mycol = scalarMap.to_rgba(float(spw))
			ax.plot(subtimes,unflaggedplot[:,:,corr],'.',markersize=6,alpha=1.0,color=mycol)
			if showflags:
				ax.plot(subtimes,flaggedplot[:,:,corr],'.',markersize=6,alpha=1.0,color='magenta')
			if gotms:
				ax.set_title(doplot+':'+str(ant)+':'+antnames[ant]+':spw'+str(spw))
			else:
				ax.set_title(str(ant)+':'+str(spw))
			if doplot in ['amp','real','imag']:
				ax.axhspan(0.9,0.9,0,1,color='black',alpha=0.2)
				ax.axhspan(1.0,1.0,0,1,color='black',alpha=0.2)
				ax.axhspan(1.1,1.1,0,1,color='black',alpha=0.2)
			elif doplot == 'phase':
				mkr = numpy.pi/2.0
				ax.axhspan(0.0,0.0,0,1,color='black',alpha=0.2)
				ax.axhspan(mkr,mkr,0,1,color='black',alpha=0.2)
				ax.axhspan(-mkr,-mkr,0,1,color='black',alpha=0.2)
		ax.set_ylim((ymin,ymax))
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.spines['left'].set_visible(False)
		ax.get_xaxis().set_visible(False)
		ax.get_yaxis().set_visible(False)
		c += 1
fig.tight_layout()
fig.savefig(pngname)
gi('Rendered '+pngname)
t.done()
