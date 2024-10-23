# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 14:21:34 2024
@author: Will Maslanka

Script for scraping data from University of Wyoming Radiosonde Page
https://weather.uwyo.edu/upperair/sounding.html

Downloads all available radiosonde ascents between two datetimes for a number 
of user-defined stations.

Main Dependancies - pandas, siphon
conda install -c conda-forge siphon
"""
import warnings
import os
from datetime import datetime
import pandas as pd
from siphon.simplewebservice.wyoming import WyomingUpperAir

warnings.filterwarnings('ignore')

def make_directory(save_dir):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

#%% User defined variables
sonde_dir = 'save-directory-here' # replace this directory with save location

# If directory doesn't exist, make it so
make_directory(sonde_dir)

stations = ['71119', '03882'] # list of WMO station codes that you want
# List of codes here: https://www.umr-cnrm.fr/dbfastex/instruments/rsc_dat.html
# Map of sites here: https://weather.uwyo.edu/upperair/sounding.html

start_dt = datetime(2024,10,1,0,0) # start datetime
end_dt = datetime(2024,10,7,0,0) # end datetime

#%% Downloading data
MAX_ATTEMPTS = 2 
# sometimes connection to archive drops out.
# Two attempts seems to work.

for att_no in range(MAX_ATTEMPTS):
    for station in stations:
        stn_dir = f'{sonde_dir}{station}/'
        make_directory(stn_dir)
        
        start_pd = pd.to_datetime(start_dt)
        end_pd = pd.to_datetime(end_dt).round(freq='12h')
        
        while start_pd <= end_pd:
            print(f'..{start_pd}')
            try:
                df = WyomingUpperAir.request_data(start_pd, station)
                df.to_csv(f"{stn_dir}{start_pd.strftime('%y%m%d%H')}_radiosonde.csv", index=False)
            except:
                pass
                    
            start_pd += pd.DateOffset(hours=12)


