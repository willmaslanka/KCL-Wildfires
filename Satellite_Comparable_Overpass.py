# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 13:15:57 2023
@author: Will Maslanka

Script for finding satellite overpasses at a given location within 15minutes.
User Defined Variables:
    earth_loc: (lat/lon location on earth)
    period of time to look at
    view angle from satellite (90deg is overhead, 0deg is on horizon)
"""

import pandas as pd
import numpy as np
import datetime
from skyfield.api import load, wgs84, N, E

#%% User Defined Functions
def find_altitude_deg(earth_loc, satellite, time):
    difference = satellite - earth_loc
    topocentric = difference.at(time)
    alt, __, __ = topocentric.altaz()
    return alt.degrees

def df_simple(df, satellite_name, event):
    df_sat = df[df.Satellite == satellite_name]
    df_sat = df_sat[df_sat.Event == event]
    df_sat = df_sat.drop(columns = "Event")
    return df_sat

#%% User Defined Variables
earth_loc = wgs84.latlon(51.067222*N, -93.793056*E) # Location of interest
altitude_degrees = 50
wanted_satellites = ["SENTINEL-3", "TERRA", "AQUA", "SUOMI NPP", "NOAA 20"]

#%% Setting up Timescales
ts = load.timescale()  # Set up a time scale

dt0 = datetime.datetime.now(datetime.timezone.utc)
dt1 = dt0 + datetime.timedelta(days=10)

t0 = ts.utc(dt0.year, dt0.month, dt0.day)  # Start time
t1 = ts.utc(dt1.year, dt1.month, dt1.day)  # End time

#%% Retrieving Satellite Paths
url = "https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FO"
satellite_list = load.tle_file(url, reload=True)
satellites = []
for sel_sat in wanted_satellites:
    for satellite in satellite_list:
        if (sel_sat in satellite.name 
            and "SAR-X" not in satellite.name 
            and "SKYTERRA" not in satellite.name):
            satellites.append(satellite)

#%% Getting Data           
a = []
for satellite in satellites:
    times_list, event_list = satellite.find_events(
        earth_loc, t0, t1, altitude_degrees=altitude_degrees)
    el_list = find_altitude_deg(earth_loc, satellite, times_list)
    ev_list = [f"rise above {altitude_degrees} deg", "culminate", f"set below {altitude_degrees}"]
    a.append([[tt.utc_strftime("%Y-%m-%d %H:%M:%S"),
               ev_list[ee],
               satellite.name,
               ev] for tt, ee, ev in zip(times_list, event_list, el_list)])
    
df = pd.concat([pd.DataFrame(x, columns=["Date", "Event", "Satellite", "Elevation"]) for x in a])
df["Date"] = pd.to_datetime(df.Date)

#%%
df_sn3a = df_simple(df,'SENTINEL-3A','culminate')
df_sn3b = df_simple(df,'SENTINEL-3B','culminate')
df_noaa = df_simple(df,'NOAA 20','culminate')
df_snpp = df_simple(df,'SUOMI NPP','culminate')
df_tera = df_simple(df,'TERRA','culminate')
df_aqua = df_simple(df,'AQUA','culminate')

df_sat_list = [df_sn3a, df_sn3b, df_noaa, df_snpp, df_tera, df_aqua]
df_check = np.zeros((6,6))

df_a_cnt = 0
for df_a in df_sat_list:
    df_a_cnt += 1
    df_b_cnt = 0
    for df_b in df_sat_list:
        df_b_cnt += 1
        
        if df_a.Satellite.iloc[0] != df_b.Satellite.iloc[0] and df_check[df_a_cnt-1, df_b_cnt-1] == 0:
            df_check[df_a_cnt-1, df_b_cnt-1] += 1
            df_check[df_b_cnt-1, df_a_cnt-1] += 1
            
            df_a_10 = df_a.Date.dt.round(freq="15T")
            df_b_10 = df_b.Date.dt.round(freq="15T")
            merged = pd.merge(df_a_10, df_b_10, on="Date", how="inner")
            
            match_ups = []
            for i, row in merged.iterrows():
                isa = np.argmin(np.abs(df_a.Date - row.Date))
                isb = np.argmin(np.abs(df_b.Date - row.Date))
                match_ups.append([df_a.iloc[isa].Date, df_b.iloc[isb].Date, abs(df_a.iloc[isa]["Date"] - df_b.iloc[isb]["Date"])])
            
            delta_df = pd.DataFrame(match_ups, columns=[f"{df_a.Satellite.iloc[0]}", f"{df_b.Satellite.iloc[0]}", "Delta"])
            
            if not delta_df.empty:
                print(f"\n{df_a.Satellite.iloc[0]} vs {df_b.Satellite.iloc[0]}")
                print('Overpasses available')
                print(delta_df)
            