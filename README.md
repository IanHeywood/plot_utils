# plot_utils

A set of command line utilities for making various plots from Measurement Sets, images and calibration tables of various flavours.

# 1) plot_ants.py

Plot antenna positions from a MS. Converts lat,lon positions to UTM eastings and northings in metres, relative to a specified array centre. Only really works for MeerKAT at the moment. Requires https://pypi.python.org/pypi/utm

# 2) plot_bandpass.py

Plot CASA bandpass tables (still needs polishing)

# 3) plot_cubical.py

Plot the contents of a CubiCal gain table. Plots one panel per antenna. Error table is overlaid on gain plot as a shaded region, and multiple directions should be automatically colour-coded.

# 4) plot_gaintab.py

Plot time-dependent gain solutions produced by CASA's gaincal task. All (or a selection of) antennas are plotted on two panels, either amp and phase, or real and imag.

# 5) plot_model.py

Plots the model from a wsclean run, overlaid on the MFS image. For visually checking the model for completeness and artefacts during self-calibration.

# 6) plot_vla_gaintab.py

This is geared towards VLA calibration but should work for any gain table, although those with only a single SPW might not make economical use of the plot space. The resulting plots have antennas in rows and SPWs in columns. The goal is to provide a global overview of the solutions to indentify misbehaving antennas or SPWs, and gauge solution coherence. This comes at the expense of any axis labels whatsoever. 

# Output gallery:

<p align="center">
  
![](https://i.imgur.com/eDzd6kK.jpg)

![](https://i.imgur.com/7u8Jox3.jpg)

![](http://i.imgur.com/plF2K6w.jpg)

</p>
