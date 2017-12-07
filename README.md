# plot_gaintab

Plot time-dependent gain solutions produced by CASA's gaincal task. All (or a selection of) antennas are plotted on two panels, either amp and phase, or real and imag.

```
Usage: plot_gaintab.py [options] tablename

Options:
  -h, --help            show this help message and exit
  -f FIELD, --field=FIELD
                        Field ID to plot (default = 0)
  -d DOPLOT, --doplot=DOPLOT
                        Plot complex values as amp and phase (ap) or real and
                        imag (ri) (default = ap)
  -a PLOTANTS, --ant=PLOTANTS
                        Plot only this antenna, or comma-separated list of
                        antennas
  -c CORR, --corr=CORR  Correlation index to plot (usually just 0 or 1,
                        default = 0)
  --t0=T0               Minimum time to plot (default = full range)
  --t1=T1               Maximum time to plot (default = full range)
  --yu0=YU0             Minimum y-value to plot for upper panel (default =
                        full range)
  --yu1=YU1             Maximum y-value to plot for upper panel (default =
                        full range)
  --yl0=YL0             Minimum y-value to plot for lower panel (default =
                        full range)
  --yl1=YL1             Maximum y-value to plot for lower panel (default =
                        full range)
  --cmap=MYCMAP         Matplotlib colour map to use for antennas (default =
                        coolwarm)
  --size=MYSIZE         Font size for figure labels (default = 20)
  -p PNGNAME, --plotname=PNGNAME
                        Output PNG name (default = something sensible)
```

![](https://i.imgur.com/eDzd6kK.jpg)

# plot_vla_gaintab

This is geared towards VLA calibration but should work for any gain table, although those with only a single SPW might not make economical use of the plot space. Uses pyrap.

The resulting plots have antennas in rows and SPWs in columns. The goal is to provide a global overview of the solutions to indentify misbehaving antennas or SPWs, and gauge solution coherence. This comes at the expense of any axis labels whatsoever. 

```
Usage: plot_casa_gaintab.py [options] gaintable

Options:
  -h, --help            show this help message and exit
  -p DOPLOT, --plot=DOPLOT
                        Set to amp, phase, real or imag (default = amp)
  -c CORR, --corr=CORR  Correlation product to plot (default = 0)
  --ms=MSNAME           Corresponding Measurement Set from which to read
                        antenna info (not essential)
  --showflags           Show flagged points in magenta (default = True)
  --y0=Y0               Minimum value of y-axis
  --y1=Y1               Maximum value of y-axis
  --pngname=PNGNAME     Name of output PNG
```

![](http://i.imgur.com/plF2K6w.jpg)
