#!/usr/bin/python

# Plot killMS gain tables
# ian.heywood@physics.ox.ac.uk


import matplotlib
matplotlib.use('Agg')
from optparse import OptionParser
import numpy
import sys
import pylab
from astropy.time import Time
from astropy.coordinates import ICRS
from astropy import units as u
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
#   ax.axes.get_yaxis().set_visible(False)


def set_fontsize(fig,fontsize):
    def match(artist):
        return artist.__module__ == 'matplotlib.text'
    for textobj in fig.findobj(match=match):
        textobj.set_fontsize(fontsize)

def time_to_index(treq,tlist):
    idx = (numpy.abs(tlist - treq)).argmin()
    return idx,tlist[idx]



def main():


# ---------------------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------------------


    parser = OptionParser(usage='%prog [options] tablename')
    parser.add_option('--list',dest='listtab',help='List table contents (default = False)',action='store_true',default=False)
    #parser.add_option('-f','--field',dest='field',help='Field ID to plot (default = 0)',default=0)
    parser.add_option('--plot',dest='doplot',help='Plot (a)mplitude, (p)hase, (r)eal or (i)maginary component of complex gains (default = a)',default='a')
    #parser.add_option('--jones',dest='jones',help='Jones matrix to plot (default = G:gain)',default='G:gain')
    parser.add_option('--ant',dest='plotants',help='Plot only this antenna, or comma-separated list of antennas',default=-1)
    parser.add_option('--dir',dest='plotdir',help='Direction (cluster) index to plot (default = all)',default=-1)
    parser.add_option('--freq',dest='plotfreq',help='Frequency index to plot (leave unset to average all)',default=-1)
    parser.add_option('--corr1',dest='corr1',help='First correlation index to plot (default = 0)',default=0)
    parser.add_option('--corr2',dest='corr2',help='Second correlation index to plot (default = 0)',default=0)
    parser.add_option('--tmin',dest='tmin',help='Minimum time to plot for all panels (default = full range)',default=-1)
    parser.add_option('--tmax',dest='tmax',help='Maximum time to plot for all panels (default = full range)',default=-1)
    parser.add_option('--ymin',dest='ymin',help='Minimum y-value to plot for all panels (default = full range per panel)',default=-1)
    parser.add_option('--ymax',dest='ymax',help='Maximum y-value to plot for all panels (default = full range per panel)',default=-1)
    parser.add_option('--unwrap',dest='unwrap',help='Unwrap phases for phase plots (default = True)',action='store_false',default=True)
    parser.add_option('--cmap',dest='mycmap',help='Matplotlib colour map to use for antennas (default = jet)',default='jet')
    parser.add_option('--size',dest='mysize',help='Font size for figure labels (default = 32)',default=32)
    parser.add_option('--ncols',dest='ncols',help='Number of columns on plot figure (default = 8)',default=8)
    parser.add_option('--msize',dest='msize',help='Marker size on plots (default = 9)',default=9)
    parser.add_option('--lwidth',dest='lwidth',help='Line width on plots (default = 2)',default=2)
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
    msize = options.msize 
    lwidth = options.lwidth 


    # ---------------------------------------------------------------------------------------
    # Some light error trapping
    # ---------------------------------------------------------------------------------------


    if len(args) != 1:
        print('Please specify a killMS table to plot.')
        sys.exit()
    else:
        gaintab = args[0]


    if doplot not in ['a','p','r','i']:
        print('Plot selection must be one of [a,p,r,i]')
        sys.exit()


    # ---------------------------------------------------------------------------------------


    tab = numpy.load(gaintab)

    clusters = tab['ClusterCat']
    freqs = tab['FreqDomains']
    ant_names = tab['StationNames']
    t0 = tab['Sols']['t0']
    t1 = tab['Sols']['t1']
    dt = t1 - t0
    t = t0 + (dt/2)
    t = t - t[0] # Everthing relative to start of observation

    gains = tab['Sols']['G']

    gshape = gains.shape
    nt = gshape[0]
    nfreq = gshape[1]
    nant = gshape[2]
    ndir = gshape[3]
    corr_i = gshape[4]
    corr_j = gshape[5]


    # ---------------------------------------------------------------------------------------
    # List the contents
    # ---------------------------------------------------------------------------------------


    if listtab:
        print('')
        print(gaintab)
        print('')

        print('    Number of antennas   '+str(nant))
        print('')

        start_time = t0[0]/86400
        end_time = t1[-1]/86400
        start_time_iso = Time(start_time,format='mjd').iso
        end_time_iso = Time(end_time,format='mjd').iso
        print('    Number of intervals  '+str(nt))
        print('    Median interval      '+str(numpy.median(dt))+' s')
        print('    Start time           '+str(t0[0])+' ('+start_time_iso+')')
        print('    End time             '+str(t1[-1])+' ('+end_time_iso+')')
        print('')
        
        print('    Frequency chunks    ',)
        for ff in freqs:
            f0 = str(round(ff[0]/1e6,2))
            f1 = str(round(ff[1]/1e6,2))
            print(f0+' - '+f1+' MHz')
            print('                        ',)
        print('')

        print('    Correlations        ',str(corr_i)+'x'+str(corr_j))
        print('')

        print('    Directions          ',)
        for i in range(0,len(clusters)):
            cl = clusters[i]
            cldir = ICRS(cl[1]*u.rad, cl[2]*u.rad)
            ra = cldir.ra.value
            dec = cldir.dec.value
            ra_hms = cldir.ra.to_string(u.hour)
            dec_dms = cldir.dec.to_string(u.deg)
            flux = str(round(cl[5],3))+' Jy'
            print('%s   %6.6f  %6.6f  %16s  %16s  %10s' % (str(i),ra,dec,ra_hms,dec_dms,flux))
            print('                        ',)
        print('')

        sys.exit()


    # ---------------------------------------------------------------------------------------
    #  Solution selection (frequency not yet implemented)
    # ---------------------------------------------------------------------------------------


    # Antenna selection


    nants = len(ant_names)
    if plotants == -1:
        ant_list = numpy.arange(0,nants)
        ant_label = 'allants'
    else:
        ant_list = []
        ant_label = 'ants'+plotants.replace(',','-')
        for p in plotants.split(','):
            p = int(p)
            if p < nants:
                ant_list.append(p)
            else:
                print('Ignoring out-of-range antenna '+str(p))


    # Direction selection


    if plotdir == -1:
        plot_dirs = numpy.arange(0,ndir)
        dir_label = 'alldirs'
    else:
        plot_dirs = []
        dir_label = 'dirs'+plotdir.replace(',','-')
        for d in plotdir.split(','):
            d = int(d)
            if d in numpy.arange(0,ndir):
                plot_dirs.append(d)
            else:
                print('Ignoring out-of-range direction '+str(d))

    cNorm = colors.Normalize(vmin=0,vmax=ndir-1)
    mymap = cm = pylab.get_cmap(mycmap)
    scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=mymap)


    # Time selection


    if tmin != -1:
        t_idx0 = time_to_index(tmin,t)[0]
    else:
        t_idx0 = 0

    if tmax != -1:
        t_idx1 = time_to_index(tmax,t)[0]
    else:
        t_idx1 = len(t)

    if t_idx1 < t_idx0:
        print('Please check your time ranges')
        sys.exit()


    # ---------------------------------------------------------------------------------------
    # Make the plots
    # ---------------------------------------------------------------------------------------


    # Generate PNG name if one isn't provided


    if pngname == '':
        pngname = 'plot_'+gaintab.replace('/','-')+'_corr'+str(corr1)+'-'+str(corr2)+'_'+doplot+'_'+ant_label+'_'+dir_label+'.png'


    figx = int(ncols*15)
    nrows = int(numpy.ceil(float(len(ant_list))/float(ncols)))
    figy = int(nrows*10)
    fig = pylab.figure(figsize=(figx,figy))

    pltcount = 1

    for ant in ant_list:

        ax1 = fig.add_subplot(nrows,ncols,pltcount,facecolor='#EEEEEE')

        
        # Label subplot with antenna / corr info
        

        plotlabel = str(ant)+':'+ant_names[ant]+' / '+str(corr1)+':'+str(corr2)
        ax1.text(0.5,1.02,plotlabel,size=mysize,horizontalalignment='center',color='black',transform=ax1.transAxes)

        
        # Initialise axis limits
        
        
        x0 = y0 = 1e20
        x1 = y1 =-1e20


        # Setup 

        
        setup_plot(ax1)


        # Loop over directions


        for mydir in plot_dirs:


            # Set colour for direction


            y1col = scalarMap.to_rgba(float(mydir))


            # Get the selected gains and time

            if plotfreq == -1:
                gains_freqavg = numpy.mean(gains,axis=1)
                g0 = (gains_freqavg[t_idx0:t_idx1,ant,mydir,corr1,corr2])
            else:
                g0 = (gains[t_idx0:t_idx1,plotfreq,ant,mydir,corr1,corr2])
            
            time = t[t_idx0:t_idx1] 


            # Amplitudes


            if doplot == 'a':
                yy = numpy.abs(g0)
                ax1.plot(time,yy,'.-',markersize=msize,lw=lwidth,alpha=1.0,zorder=150,color=y1col)


            # Phases


            elif doplot == 'p':
                yy = numpy.angle(g0)
                yy = numpy.array(yy)
                if unwrap:
                    yy = numpy.unwrap(yy)
                ax1.plot(time,yy,'.-',markersize=msize,lw=lwidth,alpha=1.0,zorder=150,color=y1col)


            # Real


            elif doplot == 'r':
                yy = numpy.real(g0)
                ax1.plot(time,yy,'.-',markersize=msize,lw=lwidth,alpha=1.0,zorder=150,color=y1col)


            # Imaginary


            elif doplot == 'i':
                yy = numpy.imag(g0)
                ax1.plot(time,yy,'.-',markersize=msize,lw=lwidth,alpha=1.0,zorder=150,color=y1col)


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
        

            if numpy.min(yy) < y0:
                if ymin == -1:
                    y0 = numpy.min(yy)
                else:
                    y0 = ymin
        

            if numpy.max(yy) > y1:
                if ymax == -1:
                    y1 = numpy.max(yy)
                else:
                    y1 = ymax


        # Apply plot limits
        

        ax1.set_xlim((x0,x1))
        ax1.set_ylim((y0,y1))


        # Final row gets the x-labels


        if pltcount > len(ant_list) - ncols:    
            for tick in ax1.get_xticklabels():
                tick.set_rotation(90)
            ax1.set_xlabel('Time - '+str(t0[0])+' [s]')
        else:
            ax1.tick_params(labelbottom='off')    


        # Left hand column gets the y-labels


        leftplots = []
        for i in range(0,ncols):
            leftplots.append((i*(nrows-1))+(i+1))
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


if __name__ == '__main__':
    main()
