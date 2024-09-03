# -*- coding: utf-8 -*-
"""
Loading Radiosonde Data from https://weather.uwyo.edu/upperair/sounding.html
Look at the website to get the station code (might be in letters or numbers)
Times are always in UTC (also known as Zulu time), and radiosondes launch
either at 00Z or 12Z.

The website can be a bit slow/fall over every so often, so if you get told that
there is a issue/delay, just wait a minute or so and try again.

Information about Radiosondes, and how they are plotted (Tephigrams, or Skew-T 
plots) can be found here: 
    https://homepages.see.leeds.ac.uk/~chmjbm/arran/radiosondes.pdf  

The siphon package info is taken from here:
    https://unidata.github.io/siphon/latest/examples/upperair/Wyoming_Request.html

@author: k2262276
"""
import pandas as pd
from siphon.simplewebservice.wyoming import WyomingUpperAir

time_wnt = pd.to_datetime("2024-09-03 00:00")
station = 'WSE'

df = WyomingUpperAir.request_data(time_wnt, station)
