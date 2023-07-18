# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:36:35 2023
@author: Will Maslanka

Script for scraping data from University of Wyoming Radiosonde Page
https://weather.uwyo.edu/upperair/sounding.html
Plots Skew-T and Hodograph, and saves to file

Main Dependancies - numpy, matplotlib, metpy, siphon
conda install -c conda-forge siphon
conda install -c conda-forge metpy
"""

from siphon.simplewebservice.wyoming import WyomingUpperAir
from metpy.plots import SkewT, Hodograph
from metpy.units import pandas_dataframe_to_unit_arrays, units
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from datetime import datetime
import metpy.calc as mpcalc
import numpy as np
import warnings

# Preamble and defining save directory
warnings.filterwarnings('ignore')
plt.close('all')

plt_dir = 'C:/Users/k2262276/Documents/FENIX_Plot/'

#%%
# Get most recent datetime (either 00z or 12z)
dt_nowutc = datetime.utcnow()
if dt_nowutc.hour < 12 and dt_nowutc.hour != 0:
    dt_wnt = datetime(dt_nowutc.year,
                      dt_nowutc.month,
                      dt_nowutc.day,
                      0,
                      0)
else:
    dt_wnt = datetime(dt_nowutc.year,
                      dt_nowutc.month,
                      dt_nowutc.day,
                      12,
                      0) 

#%% Retrieve Data
station = 'WPL' # Pickle Lake

# Read remote sounding data based on time (dt) and station
df = WyomingUpperAir.request_data(dt_wnt, station)
sounding = pandas_dataframe_to_unit_arrays(df)
#%%
# Setting Values for Mixing Ratio
w = np.array([0.01, 0.005, 0.003, 
              0.002, 0.001, 0.0006, 
              0.0002, 0.0001, 0.00005])[:, None] * units('g/g')
p = units.hPa * np.linspace(1050, 400, 9)

# Plotting Skew-T
fig = plt.figure(figsize=(12, 12))
# Grid for plots
gs = gridspec.GridSpec(3, 3)
skew = SkewT(fig, rotation=45, subplot=gs[:, :2])

# Plot the sounding using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
skew.plot(sounding['pressure'], sounding['temperature'], 'tab:red')
skew.plot(sounding['pressure'], sounding['dewpoint'], 'tab:green')

# Mask barbs to be below 100 hPa only
mask = sounding['pressure'] >= 110 * units.hPa
skew.plot_barbs(sounding['pressure'][mask][::3],
                sounding['u_wind'][mask][::3],
                sounding['v_wind'][mask][::3])

# Add the relevant special lines
skew.plot_dry_adiabats(t0=np.arange(223, 533, 10)*units.K,
                       alpha=0.25, color='orangered')
skew.plot_moist_adiabats(t0=np.arange(233, 400, 5)*units.K,
                         alpha=0.25, color='tab:green')
skew.plot_mixing_lines(mixing_ratio=w, pressure=p,
                       linestyle='dotted', color='tab:blue')

# Good bounds for aspect ratio
skew.ax.set_xlim(-30, 30)
skew.ax.set_ylim(1000, 100)

# Add Mixing Ration values
for val in w.flatten()[::1]:
    top_p = p[-1]
    dewpt = mpcalc.dewpoint(mpcalc.vapor_pressure(top_p, val))
    skew.ax.text(dewpt, top_p, str(val.to('g/kg').m),
                 horizontalalignment='center')

# Add Title to Skew-T    
plt.title('{} Sounding'.format(station), loc='left')
plt.title(f'Valid Time: {dt_wnt.strftime("%Y-%m-%d %H")}Z', loc='right')
    
# Create a hodograph
agl = sounding['height'] - sounding['height'][0]
mask = agl <= 10 * units.km
intervals = np.array([0, 1, 3, 5, 8]) * units.km
colors = ['tab:red', 'tab:green', 'tab:blue', 'tab:olive']
ax = fig.add_subplot(gs[0, -1])
h = Hodograph(ax, component_range=30.)
h.add_grid(increment=10)
l = h.plot_colormapped(sounding['u_wind'][mask],
                       sounding['v_wind'][mask],
                       agl[mask], intervals=intervals, colors=colors)
cb=plt.colorbar(l)
cb.set_label('mAGL [m]')
ax.set_xlabel('Knots')
plt.title('{} Hodograph'.format(station), loc='left')

plt.savefig(f"{plt_dir}{dt_wnt.strftime('%y%m%d%H')}Z_{station}_skewt.png",
            bbox_inches='tight', dpi=300, pad_inches = 0)
plt.close()