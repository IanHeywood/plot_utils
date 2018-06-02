#!/usr/bin/env python

# Usage: plot_ants.py <msname>
#
# Takes antenna coordinates using latitude and longitude
# and converts these to UTM eastings and northings relative
# to an array centre coordinate in lat and lon.
#
# Only really works for MeerKAT at the moment, but could add
# other antenna coordinate functions as per the meerkat_antennas()
# function, providing the 'pad' names match the NAME column of the
# antenna table in the MS.


import utm
import numpy
import sys
import pylab
from pyrap.tables import table


mysize = 20


def meerkat_antennas():
    # http://public.ska.ac.za/meerkat
    # MeerKAT64v36.wgs84.64x4.xlsx
    # Nominal array center position: 21.44389 -30.71317
    # longitude,latitude,pad
    lon0 = 21.44389
    lat0 = -30.71317
    mk_coords = [(21.44380306, -30.71292524, 'M000'),
        (21.44390085, -30.71260539, 'M001'),
        (21.44355442, -30.71307843, 'M002'),
        (21.44319474, -30.71288049, 'M003'),
        (21.44259918, -30.71333651, 'M004'),
        (21.44282392, -30.71360891, 'M005'),
        (21.44369852, -30.71372031, 'M006'),
        (21.4429543, -30.71468778, 'M007'),
        (21.44291274, -30.71588095, 'M008'),
        (21.44422681, -30.71440232, 'M009'),
        (21.444809, -30.71567238, 'M010'),
        (21.44476593, -30.7142314, 'M011'),
        (21.44535056, -30.71437675, 'M012'),
        (21.44636072, -30.71460427, 'M013'),
        (21.44681895, -30.71363266, 'M014'),
        (21.44608845, -30.71303153, 'M015'),
        (21.44689743, -30.71273217, 'M016'),
        (21.44597304, -30.71206826, 'M017'),
        (21.44499267, -30.71327305, 'M018'),
        (21.44567218, -30.71362797, 'M019'),
        (21.44490234, -30.71375803, 'M020'),
        (21.44079994, -30.71400721, 'M021'),
        (21.44052529, -30.71233809, 'M022'),
        (21.43999592, -30.71105094, 'M023'),
        (21.44022471, -30.70970248, 'M024'),
        (21.44198987, -30.70902091, 'M025'),
        (21.44285617, -30.71090172, 'M026'),
        (21.44431225, -30.71126445, 'M027'),
        (21.44335543, -30.71184184, 'M028'),
        (21.44296272, -30.71217469, 'M029'),
        (21.44567669, -30.7100276, 'M030'),
        (21.44646253, -30.71021016, 'M031'),
        (21.44870382, -30.7094729, 'M032'),
        (21.44995027, -30.70326373, 'M033'),
        (21.44762369, -30.71131072, 'M034'),
        (21.44792048, -30.71268726, 'M035'),
        (21.44794236, -30.71367768, 'M036'),
        (21.447859, -30.71519796, 'M037'),
        (21.44611607, -30.71618829, 'M038'),
        (21.44653843, -30.71639649, 'M039'),
        (21.44360936, -30.71747896, 'M040'),
        (21.440888, -30.71702307, 'M041'),
        (21.44011369, -30.71520674, 'M042'),
        (21.43731453, -30.71221265, 'M043'),
        (21.43453602, -30.70564016, 'M044'),
        (21.42475874, -30.70864938, 'M045'),
        (21.42857562, -30.69525542, 'M046'),
        (21.43785295, -30.71572093, 'M047'),
        (21.4146123, -30.68682094, 'M048'),
        (21.40625266, -30.70711438, 'M049'),
        (21.42246555, -30.71866252, 'M050'),
        (21.43501372, -30.7179944, 'M051'),
        (21.43769653, -30.72141485, 'M052'),
        (21.44398726, -30.72281975, 'M053'),
        (21.45299144, -30.71556328, 'M054'),
        (21.45643282, -30.7101845, 'M055'),
        (21.46057174, -30.70684555, 'M056'),
        (21.44696358, -30.68165598, 'M057'),
        (21.4731677, -30.68682094, 'M058'),
        (21.48236368, -30.7042055, 'M059'),
        (21.47958867, -30.72764899, 'M060'),
        (21.44371766, -30.73201279, 'M061'),
        (21.42884892, -30.73363457, 'M062'),
        (21.40819133, -30.72764899, 'M063')]
    return mk_coords,lon0,lat0


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


def plot_ants(msname):

    msname = msname.rstrip('/')

    ant_tab = table(msname+'/ANTENNA',ack=False)
    ms_ants = ant_tab.getcol('NAME')
    ms_ants = [x.upper() for x in ms_ants]
    ant_tab.close()

    ants,lon0,lat0 = meerkat_antennas()

    e0,n0,zone_num,zone_Letter = utm.from_latlon(lat0,lon0)

    ee = []
    nn = []
    labs = []

    ee_abs = []
    nn_abs = []
    labs_abs = []

    for ant in ants:
        easting,northing,zone_num,zone_letter = utm.from_latlon(ant[1],ant[0])
        if ant[2].upper() in ms_ants:
            ee.append(easting-e0)
            nn.append(northing-n0)
            labs.append(ant[2])
        else:
            ee_abs.append(easting-e0)
            nn_abs.append(northing-n0)
            labs_abs.append(ant[2])

    fig = pylab.figure(figsize=(18,18))
    ax = fig.add_subplot(1,1,1,facecolor='#EEEEEE')
    setup_plot(ax)
    ax.plot(ee,nn,'o',markersize=12,color='red',zorder=200)
    ax.plot(ee_abs,nn_abs,'o',markersize=12,color='grey',alpha=0.4,zorder=100)
    fig.suptitle(msname)
    set_fontsize(fig,mysize)

    for i in range(0,len(labs)):
        if ee[i] < 0.0:
            ax.text(ee[i]-70,nn[i],labs[i],size='x-small',verticalalignment='center',horizontalalignment='right',alpha=0.9,zorder=150)
        elif ee[i] > 0.0:
            ax.text(ee[i]+70,nn[i],labs[i],size='x-small',verticalalignment='center',horizontalalignment='left',alpha=0.9,zorder=150)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')

    pngname = 'plot_ants_'+msname.split('/')[-1]+'.png'

    fig.savefig(pngname,bbox_inches='tight')

    print 'Rendered: ',pngname

msname = sys.argv[1]
plot_ants(msname)
