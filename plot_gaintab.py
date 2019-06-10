#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from pyrap.tables import table
from optparse import OptionParser
import matplotlib.colors as colors
import matplotlib.cm as cmx
import glob
import numpy
import pylab
import sys


def setup_plot(ax):
	ax.grid(b=True,which='minor',color='white',linestyle='-',lw=2)
	ax.grid(b=True,which='major',color='white',linestyle='-',lw=2)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='x',which='both',bottom='off',top='off')
	ax.tick_params(axis='y',which='both',left='off',right='off')


def set_fontsize(fig,fontsize):
	def match(artist):
		return artist.__module__ == 'matplotlib.text'
	for textobj in fig.findobj(match=match):
		textobj.set_fontsize(fontsize)


parser = OptionParser(usage='%prog [options] tablename')
parser.add_option('-f','--field',dest='field',help='Field ID to plot (default = 0)',default=0)
parser.add_option('-d','--doplot',dest='doplot',help='Plot complex values as amp and phase (ap) or real and imag (ri) (default = ap)',default='ap')
parser.add_option('-a','--ant',dest='plotants',help='Plot only this antenna, or comma-separated list of antennas',default=[-1])
parser.add_option('-c','--corr',dest='corr',help='Correlation index to plot (usually just 0 or 1, default = 0)',default=0)
parser.add_option('--t0',dest='t0',help='Minimum time to plot (default = full range)',default=-1)
parser.add_option('--t1',dest='t1',help='Maximum time to plot (default = full range)',default=-1)
parser.add_option('--yu0',dest='yu0',help='Minimum y-value to plot for upper panel (default = full range)',default=-1)
parser.add_option('--yu1',dest='yu1',help='Maximum y-value to plot for upper panel (default = full range)',default=-1)
parser.add_option('--yl0',dest='yl0',help='Minimum y-value to plot for lower panel (default = full range)',default=-1)
parser.add_option('--yl1',dest='yl1',help='Maximum y-value to plot for lower panel (default = full range)',default=-1)
parser.add_option('--cmap',dest='mycmap',help='Matplotlib colour map to use for antennas (default = coolwarm)',default='coolwarm')
parser.add_option('--size',dest='mysize',help='Font size for figure labels (default = 20)',default=20)
parser.add_option('--ms',dest='myms',help='Measurement Set to consult for proper antenna names',default='')
parser.add_option('-p','--plotname',dest='pngname',help='Output PNG name (default = something sensible)',default='')
(options,args) = parser.parse_args()


field = int(options.field)
doplot = options.doplot
plotants = options.plotants
corr = int(options.corr)
t0 = float(options.t0)
t1 = float(options.t1)
yu0 = float(options.yu0)
yu1 = float(options.yu1)
yl0 = float(options.yl0)
yl1 = float(options.yl1)
mycmap = options.mycmap
mysize = int(options.mysize)
myms = options.myms
pngname = options.pngname


if len(args) != 1:
	print 'Please specify a gain table to plot.'
	sys.exit(-1)
else:
	mytab = args[0].rstrip('/')


if pngname == '':
	pngname = 'plot_'+mytab.split('/')[-1]+'_corr'+str(corr)+'_'+doplot+'_field'+str(field)+'.png'


if doplot not in ['ap','ri']:
	print 'Plot selection must be either ap (amp and phase) or ri (real and imag)'
	sys.exit(-1)


tt = table(mytab,ack=False)
ants = numpy.unique(tt.getcol('ANTENNA1'))
fields = numpy.unique(tt.getcol('FIELD_ID'))

cNorm = colors.Normalize(vmin=0,vmax=len(ants)-1)
mymap = cm = pylab.get_cmap(mycmap)
scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=mymap)


if int(field) not in fields.tolist():
	print 'Field ID '+str(field)+' not found'
	sys.exit(-1)


if plotants[0] != -1:
	plotants = plotants.split(',')
	for ant in plotants:
		if int(ant) not in ants:
			plotants.remove(ant)
			print 'Requested antenna ID '+str(ant)+' not found'
	if len(plotants) == 0:
		print 'No valid antennas have been requested'
		sys.exit(-1)
	else:
		plotants = numpy.array(plotants,dtype=int)
else:
	plotants = ants


if myms != '':
	anttab = table(myms.rstrip('/')+'/ANTENNA')
	antnames = anttab.getcol('NAME')
	anttab.done()
else:
	antnames = ''


fig = pylab.figure(figsize=(24,18))
ax1 = fig.add_subplot(211,facecolor='#EEEEEE')
ax2 = fig.add_subplot(212,facecolor='#EEEEEE')


setup_plot(ax1)
setup_plot(ax2)


xmin = 1e20
xmax = -1e20
ylmin = 1e20
ylmax = -1e20
yumin = 1e20
yumax = -1e20


for ant in plotants:
	y1col = scalarMap.to_rgba(float(ant))
	y2col = scalarMap.to_rgba(float(ant))

	mytaql = 'ANTENNA1=='+str(ant)+' && '
	mytaql+= 'FIELD_ID=='+str(field)

	subtab = tt.query(query=mytaql)
	cparam = subtab.getcol('CPARAM')
	flagcol = subtab.getcol('FLAG')
	times = subtab.getcol('TIME')
	times = times - times[0]
	masked_data = numpy.ma.array(data=cparam,mask=flagcol)

	if doplot == 'ap':
		y1 = numpy.abs(masked_data)[:,:,corr]
		y2 = numpy.angle(masked_data)
		y2 = numpy.array(y2[:,:,corr])
		y2 = numpy.unwrap(y2[:,0])
		ax1.plot(times,y1,'o',markersize=12,alpha=1.0,zorder=100,color=y1col)
		ax1.plot(times,y1,'-',lw=2,alpha=0.4,zorder=100,color=y1col)
		ax2.plot(times,y2,'o',markersize=12,alpha=1.0,zorder=100,color=y2col)
		ax2.plot(times,y2,'-',lw=2,alpha=0.4,zorder=100,color=y2col)
	elif doplot == 'ri':
		y1 = numpy.real(masked_data)[:,:,corr]
		y2 = numpy.imag(masked_data)[:,:,corr]
		ax1.plot(times,y1,'o',markersize=12,alpha=1.0,zorder=100,color=y1col)
		ax1.plot(times,y1,'-',lw=2,alpha=0.4,zorder=100,color=y1col)
		ax2.plot(times,y2,'o',markersize=12,alpha=0.8,zorder=100,color=y2col)
		ax2.plot(times,y2,'-',lw=2,alpha=0.4,zorder=100,color=y2col)

	subtab.close()

	dx = 1.0/float(len(ants)-1)
	if antnames == '':
		antlabel = str(ant)
	else:
		antlabel = antnames[ant]
	ax1.text(float(ant)*dx,1.05,antlabel,size='large',horizontalalignment='center',color=y1col,transform=ax1.transAxes,weight='heavy',rotation=90)

	if numpy.min(times) < xmin:
		xmin = numpy.min(times)
	if numpy.max(times) > xmax:
		xmax = numpy.max(times)
	if numpy.min(y1) < yumin:
		yumin = numpy.min(y1)
	if numpy.max(y1) > yumax:
		yumax = numpy.max(y1)
	if numpy.min(y2) < ylmin:
		ylmin = numpy.min(y2)
	if numpy.max(y2) > ylmax:
		ylmax = numpy.max(y2)


xmin = xmin-400
xmax = xmax+400
if yumin < 0.0:
	yumin = -1*(1.1*numpy.abs(yumin))
else:
	yumin = yumin*0.9
yumax = yumax*1.1
if ylmin < 0.0:
	ylmin = -1*(1.1*numpy.abs(ylmin))
else:
	ylmin = ylmin*0.9
ylmax = ylmax*1.1

if t0 != -1:
	xmin = float(t0)
if t1 != -1:
	xmax = float(t1)
if yl0 != -1:
	ylmin = yl0
if yl1 != -1:
	ylmax = yl1
if yu0 != -1:
	yumin = yu0
if yu1 != -1:
	yumax = yu1


ax1.set_xlim((xmin,xmax))
ax2.set_xlim((xmin,xmax))
ax1.set_ylim((yumin,yumax))
ax2.set_ylim((ylmin,ylmax))


if doplot == 'ap':
	ax1.set_ylabel('Amplitude')
	ax2.set_ylabel('Unwrapped phase [rad]')
elif doplot == 'ri':
	ax1.set_ylabel('Real')
	ax2.set_ylabel('Imaginary')
ax1.set_xlabel('Time [s]')
ax2.set_xlabel('Time [s]')


fig.suptitle(pngname)
set_fontsize(fig,mysize)
fig.savefig(pngname,bbox_inches='tight')

print 'Rendered: '+pngname
