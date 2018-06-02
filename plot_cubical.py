#!/usr/bin/env python
# ianh@astro.ox.ac.uk


import matplotlib
matplotlib.use('Agg')
from cubical import param_db
from optparse import OptionParser
import numpy
import sys
import pylab
import matplotlib.colors as colors
import matplotlib.cm as cmx


# ---------------------------------------------------------------------------------------
# Function definitions
# ---------------------------------------------------------------------------------------


def setup_plot(ax):
	ax.grid(b=True,which='minor',color='white',linestyle='-',lw=2)
	ax.grid(b=True,which='major',color='white',linestyle='-',lw=2)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='x',which='both',bottom='off',top='off')
	ax.tick_params(axis='y',which='both',left='off',right='off')


def dead_plot(ax):
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.tick_params(axis='x',which='both',bottom='off',top='off')
	ax.tick_params(axis='y',which='both',left='off',right='off')
	ax.axes.get_xaxis().set_visible(False)
#	ax.axes.get_yaxis().set_visible(False)


def set_fontsize(fig,fontsize):
	def match(artist):
		return artist.__module__ == 'matplotlib.text'
	for textobj in fig.findobj(match=match):
		textobj.set_fontsize(fontsize)


def gi(message):
#        print '\033[92m'+message+'\033[0m'
        print message


def ri(message):
        print '\033[91m'+message+'\033[0m'


# ---------------------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------------------


parser = OptionParser(usage='%prog [options] tablename')
parser.add_option('--list',dest='listtab',help='List parmdb contents (default = False)',action='store_true',default=False)
#parser.add_option('-f','--field',dest='field',help='Field ID to plot (default = 0)',default=0)
parser.add_option('--plot',dest='doplot',help='Plot (a)mplitude, (p)hase, (r)eal or (i)maginary component of complex gains (default = a)',default='a')
#parser.add_option('--jones',dest='jones',help='Jones matrix to plot (default = G:gain)',default='G:gain')
parser.add_option('--ant',dest='plotants',help='Plot only this antenna, or comma-separated list of antennas',default=-1)
parser.add_option('--dir',dest='plotdir',help='Direction index to plot (default = all present in selected Jones matrix)',default=-1)
parser.add_option('--freq',dest='plotfreq',help='Frequency chunk to plot (leave unset to average all)',default=-1)
parser.add_option('--corr1',dest='corr1',help='First correlation index to plot (default = 0)',default=0)
parser.add_option('--corr2',dest='corr2',help='Second correlation index to plot (default = 0)',default=0)
parser.add_option('--tmin',dest='tmin',help='Minimum time to plot for all panels (default = full range)',default=-1)
parser.add_option('--tmax',dest='tmax',help='Maximum time to plot for all panels (default = full range)',default=-1)
parser.add_option('--ymin',dest='ymin',help='Minimum y-value to plot for all panels (default = full range per panel)',default=-1)
parser.add_option('--ymax',dest='ymax',help='Maximum y-value to plot for all panels (default = full range per panel)',default=-1)
parser.add_option('--unwrap',dest='unwrap',help='Unwrap phases for phase plots (default = True)',action='store_false',default=True)
parser.add_option('--cmap',dest='mycmap',help='Matplotlib colour map to use for antennas (default = coolwarm)',default='coolwarm')
parser.add_option('--size',dest='mysize',help='Font size for figure labels (default = 32)',default=32)
parser.add_option('--ncols',dest='ncols',help='Number of columns on plot figure (default = 8)',default=8)
parser.add_option('--pngname',dest='pngname',help='Output PNG name (default = something verbose)',default='')
(options,args) = parser.parse_args()


listtab = options.listtab 
#field = int(options.field)
doplot = options.doplot
#jones = options.jones
plotants = options.plotants
plotdir = options.plotdir
plotfreq = int(options.plotfreq)
corr1 = int(options.corr1)
corr2 = int(options.corr2)
tmin = float(options.tmin)
tmax = float(options.tmax)
ymin = float(options.ymin)
ymax = float(options.ymax)
unwrap = options.unwrap
mycmap = options.mycmap
mysize = int(options.mysize)
ncols = int(options.ncols)
pngname = options.pngname


# ---------------------------------------------------------------------------------------
# Some light error trapping
# ---------------------------------------------------------------------------------------


if len(args) != 1:
	ri('Please specify a CubiCal parmdb to plot.')
	sys.exit()
else:
	gaintab = args[0].rstrip('/')


if doplot not in ['a','p','r','i']:
	ri('Plot selection must be one of [a,p,r,i]')
	sys.exit()


# ---------------------------------------------------------------------------------------


db = param_db.load(gaintab)
contents = db.names()
for item in contents:
	if item.split(':')[-1] == 'gain':
		jones = item
gain = db[jones]
gainerr = db[jones+'.err']


# ---------------------------------------------------------------------------------------
# List the contents
# ---------------------------------------------------------------------------------------


if listtab:
	gi('')
	gi('      '+gaintab)
	gi('')

	contents = db.names()
	gi('      Jones term: '+jones)
	gi('')
	interps = gain.interpolation_axes
	shape = gain.shape
	labels = gain.axis_labels
	gi('      ------------------------------')
	gi('      Axis  : Shape  :')
	gi('      ------------------------------')
	for i in range(0,len(labels)):
		if i in interps:
			info = 'Interpolable'
		elif labels[i] == 'ant':
			valid_ants = [ant for ant in range(len(gain.grid[gain.ax.ant])) 
				if gain.is_slice_valid(dir=0,ant=ant,corr1=0,corr2=0)]
			n_valid = len(valid_ants)
			info = str(n_valid)+' valid (in dir 0)'
		else:
			info = ''
		gi('      %-6s: %-6s : %-12s'%(labels[i],str(shape[i]),info))
	gi('      ------------------------------')
	gi('')
	sys.exit()


# ---------------------------------------------------------------------------------------
# 
# ---------------------------------------------------------------------------------------


# Setup which antennas to plot
ant_names = gain.grid[gain.ax.ant]
nants = len(ant_names)
valid_ants = [ant for ant in range(len(gain.grid[gain.ax.ant])) if gain.is_slice_valid(dir=0,ant=ant,corr1=0,corr2=0)]
if plotants == -1:
	ant_list = numpy.arange(0,nants)
else:
	ant_list = []
	for p in plotants.split(','):
		p = int(p)
		if p < nants:
			ant_list.append(p)
		else:
			ri('Ignoring requested out-of-range antenna index: '+str(p))


# Setup which dirs to plot, and the direction colour scale
if plotdir == -1:
	plot_dirs = gain.grid[gain.ax.dir]
else:
	plot_dirs = []
	for d in plotdir.split(','):
		d = int(d)
		if d in gain.grid[gain.ax.dir]:
			plot_dirs.append(d)
		else:
			ri('Ignoring requested out-of-range direction index: '+str(d))
ndir = len(gain.grid[gain.ax.dir])
cNorm = colors.Normalize(vmin=0,vmax=ndir-1)
mymap = cm = pylab.get_cmap(mycmap)
scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=mymap)


# Generate PNG name if one isn't provided
if pngname == '':
	if plotfreq == -1:
		fstr = 'f-avg'
	else:
		fstr = 'f'+str(plotfreq)
# 	pngname = 'plot_'+gaintab.replace('/','-')+'_corr'+str(corr1)+'-'+str(corr2)+'_'+doplot+'_field'+str(field)+'.png'
# figtitle = gaintab+' corr1='+str(corr1)+' corr2='+str(corr2)+' field='+str(field)
	pngname = 'plot_'+gaintab.replace('/','-')+'_'+fstr+'_corr'+str(corr1)+'-'+str(corr2)+'_'+doplot+'.png'


# Fig layout
figx = int(ncols*15)
nrows = int(numpy.ceil(float(len(ant_list))/float(ncols)))
figy = int(nrows*10)
fig = pylab.figure(figsize=(figx,figy))

msize = 9


pltcount = 1

for ant in ant_list:
	# Create subplot
	ax1 = fig.add_subplot(nrows,ncols,pltcount,facecolor='#EEEEEE')

	# Label subplot with antenna / corr info
	plotlabel = str(ant)+':'+ant_names[ant]+' / '+str(corr1)+':'+str(corr2)
	ax1.text(0.5,1.02,plotlabel,size=mysize,horizontalalignment='center',color='black',transform=ax1.transAxes)

	# Initialise axis limits
	x0 = y0 = 1e20
	x1 = y1 =-1e20

	# Loop over antennas
	if ant in valid_ants:

		# Setup for valid plot
		setup_plot(ax1)

		# Loop over directions
		for mydir in plot_dirs:

			# Set colour for direction
			y1col = scalarMap.to_rgba(float(mydir))


			# Get gain info from relevant slice
			g0, (t, freq) = gain.get_slice(dir=mydir,ant=ant,corr1=corr1,corr2=corr2)
			ge, (t, freq) = gainerr.get_slice(dir=mydir,ant=ant,corr1=corr1,corr2=corr2)

			# Select or average frequency chunks
			if plotfreq == -1:
				g0 = numpy.mean(g0,axis=1)
				ge = numpy.mean(ge,axis=1)
			else:
				g0 = g0[:,plotfreq]
				ge = ge[:,plotfreq]


			# Apply offset for label clarity
			time = t - t[0]


			# Amplitudes
			if doplot == 'a':
				yy = numpy.abs(g0)
				ye = numpy.abs(ge)
				ax1.plot(time,yy,'.-',markersize=msize,alpha=1.0,zorder=150,color=y1col)
				ax1.fill_between(time,y1=(yy-ye),y2=(yy+ye),color=y1col,zorder=100,alpha=0.1)

			# Phases
			elif doplot == 'p':
				yy = numpy.angle(g0)
				yy = numpy.array(yy)
				ye = numpy.angle(ge)
				ye = numpy.array(ye)
				if unwrap:
					yy = numpy.unwrap(yy)
					ye = numpy.unwrap(ye)
				ax1.plot(time,yy,'.-',markersize=msize,alpha=1.0,zorder=150,color=y1col)
				ax1.fill_between(time,y1=(yy-ye),y2=(yy+ye),color=y1col,zorder=100,alpha=0.1)

			# Real
			elif doplot == 'r':
				yy = numpy.real(g0)
				ye = numpy.real(ge)
				ax1.plot(time,yy,'.-',markersize=msize,alpha=1.0,zorder=150,color=y1col)
				ax1.fill_between(time,y1=(yy-ye),y2=(yy+ye),color=y1col,zorder=100,alpha=0.1)

			# Imaginary
			elif doplot == 'i':
				yy = numpy.imag(g0)
				ye = numpy.imag(ge)
				ax1.plot(time,yy,'.-',markersize=msize,alpha=1.0,zorder=150,color=y1col)
				ax1.fill_between(time,y1=(yy-ye),y2=(yy+ye),color=y1col,zorder=100,alpha=0.1)


			# Adjust plot limits depending on data / user prefs
			if numpy.min(time) < x0:
				if tmin == -1:
					x0 = numpy.min(time)
				else:
					x0 = tmin
			if numpy.max(time) > x1:
				if tmax == -1:
					x1 = numpy.max(time)
				else:
					x1 = tmax
			if numpy.min(yy-ye) < y0:
				if ymin == -1:
					y0 = numpy.min(yy)
				else:
					y0 = ymin
			if numpy.max(yy+ye) > y1:
				if ymax == -1:
					y1 = numpy.max(yy)
				else:
					y1 = ymax


		# Apply plot limits
		ax1.set_xlim((x0,x1))
		ax1.set_ylim((y0,y1))


		# X-labels for final row
		if pltcount > len(ant_list) - ncols:	
			for tick in ax1.get_xticklabels():
				tick.set_rotation(90)
			ax1.set_xlabel('Time - '+str(t[0])+' [s]')
		else:
			ax1.tick_params(labelbottom='off')    


	# Plot for a non-valid antenna
	else:
		ax1.plot([0,1],[0,1],'-',color='#EDD1D1',lw=50,alpha=1.0)
		ax1.plot([1,0],[0,1],'-',color='#EDD1D1',lw=50,alpha=1.0)
		dead_plot(ax1)
		ax1.tick_params(labelleft='off')
		# continue


	# Y-axis labels for plots on the left hand side
	leftplots = []
	for i in range(0,ncols):
		leftplots.append((i*(nrows))+(i+1))
	if pltcount in leftplots:
		if doplot == 'a':
			ax1.set_ylabel('Amplitude')
		elif doplot == 'p':
			ax1.set_ylabel('Phase [rad]')
		elif doplot == 'r':
			ax1.set_ylabel('Real')
		elif doplot == 'i':
			ax1.set_ylabel('Imaginary')


	set_fontsize(ax1,mysize)
	pltcount += 1


# fig.suptitle(figtitle)
fig.savefig(pngname,bbox_inches='tight')


