# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 09:30:34 2023

@author: k2262276
"""
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
import matplotlib.ticker as mticker
import numpy as np
import warnings
from datetime import date

#%%
warnings.filterwarnings('ignore')
plt.close('all')

#%%
def data_retrieve(sat, north, east, west, south, daynight, day_range):
    """
    Downloads FIRMS NRT FRP data via api.

    Parameters
    ----------
    sat : String
        Satellite Wanted (AQUA, TERRA, NOAA-20, SNPP).
    north : float
        Maximum Latitude.
    east : float
        Maximum Longitude.
    west : float
        Minimum Longitude.
    south : float
        Minimum Latitude.
    daynight : String
        Day Night Indicator.
    day_range: float
        How many days previous to check from now

    Returns
    -------
    df : Data Frame
        FRP data scrapped via API from FIRMS for the chosen Satellite.
    """
    url = "https://firms.modaps.eosdis.nasa.gov/api/area/csv/"
    api = "c901414b99ab783017e5d2f87133e777"
    
    if sat.upper() == 'SNPP':
        sat_url = "/VIIRS_SNPP_NRT/"
    elif sat.upper() == 'NOAA-20':
        sat_url = "/VIIRS_NOAA20_NRT/"
    elif sat.upper() == 'AQUA' or sat.upper() == 'TERRA':
        sat_url = "/MODIS_NRT/"

    fname = f"{url}{api}{sat_url}{np.floor(west)},{np.floor(south)},{np.ceil(east)},{np.ceil(north)}/{day_range}"

    df = pd.read_csv(fname)
    
    # If MODIS, filter by Satellite
    if sat.upper() == 'AQUA' or sat.upper() == 'TERRA':
        df = df[df['satellite'] == sat.title()]
    
    if daynight.upper() == 'D' or daynight.upper() == 'DAY':
        dayid = 'D'
    elif daynight.upper() == 'N' or daynight.upper() == 'NIGHT':
        dayid = 'N'
    
    
    df = df[(df['latitude'] >= south) & 
                     (df['latitude'] <= north) &
                     (df['longitude'] >= west) &
                     (df['longitude'] <= east) &
                     (df['daynight'] == dayid)]
    return df

def clipping_df(extent_df):
    """
    Clipping manually to remove non-Ontairo points
    NOT NEEDED IF LOOKING AT REGION NOT DEFINED BY TERRITORY EDGES

    Parameters
    ----------
    time_df : Data Frame
        Data Frame of FRP data.
    Returns
    -------
    mask_df : Data Frame
        Trimmed Data Frame.

    """
    to_remove = extent_df[
        (extent_df['latitude'] <= 45.5) & (extent_df['longitude'] <= -82.5) |
        (extent_df['latitude'] >= 46.0) & (extent_df['longitude'] >= -77.5) |
        (extent_df['latitude'] >= 47.0) & (extent_df['longitude'] >= -79.5) |
        (extent_df['latitude'] <= 47.5) & (extent_df['longitude'] <= -85.0) |
        (extent_df['latitude'] <= 48.0) & (extent_df['longitude'] <= -90.0) |
        (extent_df['latitude'] >= 55.00) & (extent_df['longitude'] <= -91.85) |
        (extent_df['latitude'] >= 54.55) & (extent_df['longitude'] <= -92.50) |
        (extent_df['latitude'] >= 55.70) & (extent_df['longitude'] <= -90.78) |
        (extent_df['latitude'] >= 56.40) & (extent_df['longitude'] <= -89.69) |
        (extent_df['latitude'] >= 56.00) & (extent_df['longitude'] <= -90.39)]

    keep_df = pd.merge(extent_df,to_remove, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)
    mask_df = datetime_df(keep_df)
    
    return mask_df

def datetime_df(full_df):
    """
    Ensures the Dataframe has acquisition date and time setup correctly (UTC).

    Parameters
    ----------
    full_df : Dataframe
        Clipped Dataframe of FRP data.

    Returns
    -------
    mask_df : Dataframe
        Dataframe with datetime column added.

    """
    mask_df = full_df.copy()
    mask_df['year'] = mask_df['acq_date'].str[:4]
    mask_df['month'] = mask_df['acq_date'].str[5:7]
    mask_df['day'] = mask_df['acq_date'].str[8:10]
    mask_df['acq_time_str'] = mask_df['acq_time'].astype(str)
    mask_df['hour'] = mask_df['acq_time_str'].str[:-2].replace('',0)
    mask_df['minute'] = mask_df['acq_time_str'].str[-2:]
    mask_df['acq_datetime'] = pd.to_datetime(mask_df[['year','month','day','hour','minute']])
    return mask_df

def unique_orbit(clip_df):
    """
    Calculates which FRP data is associated with a common orbit (done via time
                                                                 of acquisiton)

    Parameters
    ----------
    clip_df : Dataframe
        Dataframe with datetime column added.

    Returns
    -------
    rel_orb : List
        Indicies of associated orbits, using un array times.
    un : Array
        Array of unique timesteps within FRP data.

    """
    un = np.unique(clip_df.acq_datetime.values)

    rel_orb_id = 0
    rel_orb = [rel_orb_id]

    for en_a, tme in enumerate(un[1:]):
        if (tme-un[en_a]).astype('timedelta64[m]').astype(int) <= 5:
            rel_orb.append(rel_orb_id)
        else:
            rel_orb_id += 1
            rel_orb.append(rel_orb_id)
    
    return rel_orb, un

def plotting(rel_orb_df, north, east, west, south, pickle_lake, red_lake, srt, end, sat, plot_dir):
    """
    Plotting Function

    Parameters
    ----------
    rel_orb_df : TYPE
        DESCRIPTION.
    north : float
        Maximum Latitude.
    east : float
        Maximum Longitude.
    west : float
        Minimum Longitude.
    south : float
        Minimum Latitude.
    pickle_lake : List
        Longitude, Latitude of Pickle Lake.
    red_lake : List
        Longitude, Latitude of Red Lake.
    srt : datetime64
        First Datetime on the relative orbit.
    end : datetime64
        Last Datetime on the relative orbit.
    sat : String
        Satellite Wanted (AQUA, TERRA, NOAA-20, SNPP).
    plot_dir : String
        Directory for Plots to be saved.

    Returns
    -------
    None.

    """
    dist = dist_det(north, east, west, south)

    fig, ax = plt.subplots(figsize=(15,10),
                           subplot_kw={"projection": ccrs.Mercator()})
    ax.set_extent([west-dist, east+dist, south-dist, north+dist],
                  crs=ccrs.PlateCarree())
    sca = plt.scatter(rel_orb_df.longitude.values, rel_orb_df.latitude.values,
                      s=5,
                      c=np.log10(rel_orb_df.frp.values),
                      vmin = 0,
                      vmax = 3,
                      cmap="hot",
                      transform=ccrs.PlateCarree(),
                      zorder = 2)
    plt.scatter(pickle_lake[0], pickle_lake[1],
                marker='X',
                c='dodgerblue',
                transform=ccrs.PlateCarree(),
                zorder = 1)
    plt.scatter(red_lake[0], red_lake[1],
                marker='*',
                c='darkviolet',
                transform=ccrs.PlateCarree(),
                zorder = 1)
    
    cbar = plt.colorbar(sca)
    cbar.set_label(label = '$log_{10}(frp)$', 
                  size = 14, weight = 'bold')

    ax.add_feature(cf.BORDERS, edgecolor='black')
    ax.add_feature(cf.COASTLINE, edgecolor='black')
    ax.add_feature(cf.STATES, edgecolor='black')
    ax.set_facecolor("gainsboro")

    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.5,
        color="gray",
        alpha=1.0,
        linestyle="--",
        )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(np.arange(-95, -70, 2.5))
    gl.ylocator = mticker.FixedLocator(np.arange(40, 60, 2.5))
    
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    lon_formatter = LongitudeFormatter(number_format='.2f',
                                   degree_symbol='',
                                   dateline_direction_label=True)
    lat_formatter = LatitudeFormatter(number_format='.2f',
                                  degree_symbol='')
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    plt.title(f"{sat}\n{srt} --> {end} UTC")
    
    srt_acq = ACQ_datetime(srt, sat)
    end_acq = ACQ_datetime(end, sat)
    
    
    plt.savefig(f"{plot_dir}{srt_acq.split('.')[0]}_{srt_acq.split('.')[1]}_{end_acq.split('.')[1]}_{sat}_FRP.png",
                bbox_inches='tight', dpi=300, pad_inches = 0)
    plt.close()
    
def dist_det(north, east, west, south):
    """
    Calculates appropriate edge gap for plotting

    Parameters
    ----------
    north : float
        Maximum Latitude.
    east : float
        Maximum Longitude.
    west : float
        Minimum Longitude.
    south : float
        Minimum Latitude.

    Returns
    -------
    dist : float
        edge gap calculated.

    """
    dist = np.ceil(north)-north
    if np.ceil(west)-west > dist:
        dist = np.ceil(west)-west
    if south - np.floor(south) > dist:
        dist = south - np.floor(south)
    if east - np.floor(east) > dist:
        dist = east - np.floor(east)
        
    return dist

def ACQ_datetime(date_want, sat):
    """
    Calculates acquistation datetime for filename 

    Parameters
    ----------
    date_want : datetime64
        Time of observation.
    sat : String
        Satellite Wanted (AQUA, TERRA, NOAA-20, SNPP).

    Returns
    -------
    acq_datetime : String
        String for Filename.

    """
    year_wnt = str(date_want).split('-')[0]
    month_wnt = str(date_want).split('-')[1]
    day_wnt = str(date_want).split('-')[2].split('T')[0]
    hour_wnt = str(date_want).split('T')[1].split(':')[0]
    minute_wnt = str(date_want).split('T')[1].split(':')[1]

    doy = doy_finder(year_wnt, month_wnt, day_wnt)

    if sat.upper() == 'SNPP' or sat.upper() == 'NOAA-20' :
        min_floor = str((int(minute_wnt) // 6)*6)
    else:
        min_floor = str((int(minute_wnt) // 5)*5)

    if len(min_floor) != 2:
        min_floor = '0' + min_floor
        
    time_wnt = hour_wnt + min_floor
    
    acq_datetime = 'A' + str(year_wnt) + str(doy) + "." + time_wnt
    
    return acq_datetime

def doy_finder(yr, mn, dy):
    """Calculates the day of year from user input year, month,day.
    Input: yr: year (yyyy)
    Input: mn: month (mm)
    Input: dy: day (dd)

    Output: doy: Day of Year (1st Jan = 1)"""

    doy = date(int(yr), int(mn), int(dy)).timetuple().tm_yday
    return doy

#%%
if __name__ == "__main__":
    save_dir = 'C:/Users/k2262276/Documents/FENIX_SAT/'
    plot_dir = 'C:/Users/k2262276/Documents/FENIX_Plot/'
    
    north = 56.86
    east = -74.34
    west = -95.16
    south = 41.66
    day_range = 4
    daynight = 'D' # 'D'/'N'
    sat = 'TERRA' # 'AQUA' / 'TERRA' / 'SNPP' / 'NOAA-20'
    pickle_lake = [-90.18333, 51.466667]
    red_lake = [-93.793056, 51.0679222]
    
    #%%
    df = data_retrieve(sat, north, east, west, south, daynight, day_range)   
    clip_df = clipping_df(df)    
    rel_orb, un = unique_orbit(clip_df)
    
    #%%
    for ro_val in np.arange(np.max(rel_orb)+1):
        srt = np.datetime64(un[[i for i, j in enumerate(rel_orb) if j == ro_val]].min(),'m')
        end = np.datetime64(un[[i for i, j in enumerate(rel_orb) if j == ro_val]].max(),'m')
        
        rel_orb_df = clip_df[clip_df['acq_datetime'].isin(un[[i for i, j in enumerate(rel_orb) if j == ro_val]])]        
        plotting(rel_orb_df, north, east, west, south, pickle_lake, red_lake, srt, end, sat, plot_dir) 