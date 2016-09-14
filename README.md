# plot_casa_gaintab

Plot time-dependent gain solutions produced by CASA's gaincal task. This is geared towards VLA calibration but should work for any gain table, although those with only a single SPW might not make economical use of the plot space. Uses pyrap.

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
