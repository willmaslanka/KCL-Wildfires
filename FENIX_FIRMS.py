# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 09:30:34 2023
@author: Will Maslanka

Plots the FRP locations for a given point, within 250 km of the user define location.
FRP data from FIRMS
Satellites used, NOAA20, SNPP, AQUA, TERRA, LANDSAT

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
from datetime import date, datetime
from cartopy.geodesic import Geodesic

#%% Preamble
warnings.filterwarnings('ignore')
plt.close('all')

#%% Pre-Defined Functions
def tens_str(num):
    """ Turns float < 100  to 2 digit string """
    if int(num) >= 10:
        num_str = str(num)
    else:
        num_str = '0'+str(num)
    return num_str

def data_retrieve(sat, location, daynight, day_range, date_want):
    """ Downloads FIRMS NRT FRP data via api. """

    url = "https://firms.modaps.eosdis.nasa.gov/api/area/csv/"
    api = "c901414b99ab783017e5d2f87133e777"
    
    # Applying Satellite Filter
    if sat.upper() == 'SNPP':
        sat_url = "/VIIRS_SNPP_NRT/"
    elif sat.upper() == 'NOAA-20':
        sat_url = "/VIIRS_NOAA20_NRT/"
    elif sat.upper() == 'AQUA' or sat.upper() == 'TERRA':
        sat_url = "/MODIS_NRT/"
    elif sat.upper() == 'LANDSAT':
        sat_url = "/LANDSAT_NRT/"
    
    # If datetime.now is used to get time, then remove date_url from url
    # FIRMS doesn't like it for some reason
    if f"{datetime.now().year}-{tens_str(datetime.now().month)}-{tens_str(datetime.now().day)}" != f"{date_want.year}-{tens_str(date_want.month)}-{tens_str(date_want.day)}":
        date_url = f"{date_want.year}-{tens_str(date_want.month)}-{tens_str(date_want.day)}"    
    else:
        date_url = ""
    
    # Calculating extents of 250km radius circle
    cp250 = circle_calc(location, 250)
    west = cp250[:,0].min()
    east = cp250[:,0].max()
    north = cp250[:,1].max()
    south = cp250[:,1].max()
     
    # Reading in FIRMS data
    fname = f"{url}{api}{sat_url}{np.floor(west)},{np.floor(south)},{np.ceil(east)},{np.ceil(north)}/{day_range}/{date_url}"
    df = pd.read_csv(fname)
    
    # If MODIS, filter by Satellite
    if sat.upper() == 'AQUA' or sat.upper() == 'TERRA':
        df = df[df['satellite'] == sat.title()]
    
    # If LANDSAT, add column of FRP values = 10 (no FRP on LANDSAT data)
    if sat.upper() == 'LANDSAT':
        df['frp'] = 10
    
    # Get Day/Night ID
    if daynight.upper() == 'D' or daynight.upper() == 'DAY':
        dayid = 'D'
    elif daynight.upper() == 'N' or daynight.upper() == 'NIGHT':
        dayid = 'N'

    # Clipping check for 250km radius search
    df = df[(df['latitude'] >= cp250[:,1].min()) & 
                     (df['latitude'] <= cp250[:,1].max()) &
                     (df['longitude'] >= cp250[:,0].min()) &
                     (df['longitude'] <= cp250[:,0].max()) &
                     (df['daynight'] == dayid)]
    
    # Applying Datetime Function
    new_df = datetime_df(df)

    return new_df

def datetime_df(full_df):
    """ Ensures the Dataframe has acquisition date and time setup correctly (UTC). """

    mask_df = full_df.copy()
    mask_df['year'] = mask_df['acq_date'].str[:4]
    mask_df['month'] = mask_df['acq_date'].str[5:7]
    mask_df['day'] = mask_df['acq_date'].str[8:10]
    mask_df['acq_time_str'] = mask_df['acq_time'].astype(str)
    mask_df['hour'] = mask_df['acq_time_str'].str[:-2].replace('',0)
    mask_df['minute'] = mask_df['acq_time_str'].str[-2:]
    mask_df['acq_datetime'] = pd.to_datetime(mask_df[['year','month','day','hour','minute']])
    return mask_df

def circle_calc(pickle_lake, radius=200):
    """ Calculates circle with defined radius (default 200km) around a given point """
    gd = Geodesic()
    cp = gd.circle(lon=pickle_lake[0], lat=pickle_lake[1], radius = radius*1e3)
    cp = np.vstack((cp, cp[0,:]))
    return cp

def unique_orbit(clip_df):
    """ Calculates which FRP data is associated with a common orbit (done via 
    time of acquisiton) """

    # Get unique timestamps
    un = np.unique(clip_df.acq_datetime.values)

    rel_orb_id = 0
    rel_orb = [rel_orb_id]

    # If the timestamp is less than 5 minutes apart, it is part of the same 
    # relative orbit
    for en_a, tme in enumerate(un[1:]):
        if (tme-un[en_a]).astype('timedelta64[m]').astype(int) <= 5:
            rel_orb.append(rel_orb_id)
        else:
            rel_orb_id += 1
            rel_orb.append(rel_orb_id)
    
    return rel_orb, un

def plotting(rel_orb_df, pickle_lake, srt, end, sat, plot_dir):
    """ Plotting Function """

    # Calculate Circles with km radius
    cp250 = circle_calc(pickle_lake, 250)
    cp200 = circle_calc(pickle_lake, 200)
    cp150 = circle_calc(pickle_lake, 150)
    cp100 = circle_calc(pickle_lake, 100)
    cp050 = circle_calc(pickle_lake, 50)
    
    ## NOTE ##
    # Circle Calculation centre of circle to be slightly more north than input
    # co-ordinates. Best Solution is to calculate the offset, and apply it.
    ## END NOTE ##
    
    # circle calculated mid latitude
    cp_mid_lat = cp250[45,1] 
    
    # Add Offset to latitude 
    #(Calculates latitudes ~ 2.5km south of where it should be)
    cp250[:,1] += (pickle_lake[1]-cp_mid_lat)
    cp200[:,1] += (pickle_lake[1]-cp_mid_lat)
    cp150[:,1] += (pickle_lake[1]-cp_mid_lat)
    cp100[:,1] += (pickle_lake[1]-cp_mid_lat)
    cp050[:,1] += (pickle_lake[1]-cp_mid_lat)

    ## Plotting Starts here ##
    # Make Figure Axis and Extent
    fig, ax = plt.subplots(figsize=(15,10),
                           subplot_kw={"projection": ccrs.Mercator()})
    ax.set_extent([cp250[:,0].min()-0.1, # West
                   cp250[:,0].max()+0.1, # East
                   cp250[:,1].min()-0.1, # South
                   cp250[:,1].max()+0.1], # North
                  crs=ccrs.PlateCarree())
    
    # Plot FRP Locations, colourbar is log10(FRP)
    sca = plt.scatter(rel_orb_df.longitude.values,
                      rel_orb_df.latitude.values,
                      s=5,
                      c=np.log10(rel_orb_df.frp.values),
                      vmin = 0,
                      vmax = 3,
                      cmap="hot",
                      transform=ccrs.PlateCarree(),
                      zorder = 10)
    
    # Plot Location Centre (User Defined Lat/Lon)
    plt.scatter(pickle_lake[0],
                pickle_lake[1],
                marker='X',
                c='dodgerblue',
                edgecolor='black',
                transform=ccrs.PlateCarree(),
                zorder = 10)
        
    # Plotting radius circles  (250km 200km, 150km, 100km, 50km)      
    plt.plot(cp250[:,0], cp250[:,1],
             transform=ccrs.PlateCarree(),
             linestyle = 'dotted',
             color='darkslategrey',
             zorder = 0)
    
    plt.plot(cp200[:,0], cp200[:,1],
             transform=ccrs.PlateCarree(),
             linestyle = 'dotted',
             color='darkslategrey',
             zorder = 0)
    
    plt.plot(cp150[:,0], cp150[:,1],
             transform=ccrs.PlateCarree(),
             linestyle = 'dotted',
             color='darkslategrey',
             zorder = 0)
    
    plt.plot(cp100[:,0], cp100[:,1],
             transform=ccrs.PlateCarree(),
             linestyle = 'dotted',
             color='darkslategrey',
             zorder = 0)
    
    plt.plot(cp050[:,0], cp050[:,1],
             transform=ccrs.PlateCarree(),
             linestyle = 'dotted',
             color='darkslategrey',
             zorder = 0)
    
    # Plotting Bearings Lines
    # 0 - 180
    plt.plot([cp250[0,0], cp250[90,0]],
             [cp250[0,1], cp250[90,1]],
             transform=ccrs.PlateCarree(),
             color='darkslategrey',
             zorder = 0)
    
    # 30 - 210
    plt.plot([cp250[165,0], cp250[75,0]],
             [cp250[165,1], cp250[75,1]],
             transform=ccrs.PlateCarree(),
             color='darkslategrey',
             zorder = 0)
    
    # 60 - 240
    plt.plot([cp250[150,0], cp250[60,0]],
             [cp250[150,1], cp250[60,1]],
             transform=ccrs.PlateCarree(),
             color='darkslategrey',
             zorder = 0)
    
    # 90 - 270
    plt.plot([cp250[135,0], cp250[45,0]],
             [cp250[135,1], cp250[45,1]],
             transform=ccrs.PlateCarree(),
             color='darkslategrey',
             zorder = 0)
    
    # 120 - 300
    plt.plot([cp250[120,0], cp250[30,0]],
             [cp250[120,1], cp250[30,1]],
             transform=ccrs.PlateCarree(),
             color='darkslategrey',
             zorder = 0)
    
    # 150 - 330
    plt.plot([cp250[105,0], cp250[15,0]],
             [cp250[105,1], cp250[15,1]],
             transform=ccrs.PlateCarree(),
             color='darkslategrey',
             zorder = 0)
    
    # Labels for circles
    ax.text(pickle_lake[0], cp250[:,1].min(), '250km',
            transform=ccrs.PlateCarree(),
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='center')
    
    ax.text(pickle_lake[0], cp200[:,1].min(), '200km',
            transform=ccrs.PlateCarree(),
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='center')
    
    ax.text(pickle_lake[0], cp150[:,1].min(), '150km',
            transform=ccrs.PlateCarree(),
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='center')
    
    ax.text(pickle_lake[0], cp100[:,1].min(), '100km',
            transform=ccrs.PlateCarree(),
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='center')
    
    ax.text(pickle_lake[0], cp050[:,1].min(), '50km',
            transform=ccrs.PlateCarree(),
            fontsize=12,
            verticalalignment='top',
            horizontalalignment='center')
    
    # Colorbar
    cbar = plt.colorbar(sca)
    cbar.set_label(label = '$log_{10}(frp)$', 
                  size = 14, weight = 'bold')
    
    # Extra Features (incase location changes)
    ax.add_feature(cf.BORDERS, edgecolor='black')
    ax.add_feature(cf.COASTLINE, edgecolor='black')
    ax.add_feature(cf.STATES, edgecolor='black')

    # Adding Gridlines
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
    gl.xlocator = mticker.FixedLocator(np.arange(-95, -70, 1))
    gl.ylocator = mticker.FixedLocator(np.arange(40, 60, 1))
    
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    lon_formatter = LongitudeFormatter(number_format='.2f',
                                   degree_symbol='',
                                   dateline_direction_label=True)
    lat_formatter = LatitudeFormatter(number_format='.2f',
                                  degree_symbol='')
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    
    plt.title(f"{sat}\n{srt} UTC")
    srt_acq = ACQ_datetime(srt, sat)  
    plt.savefig(f"{plot_dir}{srt_acq.split('.')[0]}_{srt_acq.split('.')[1]}_{sat}_FRP.png",
                bbox_inches='tight', dpi=300, pad_inches = 0)
    plt.close()
    print(f"{plot_dir}{srt_acq.split('.')[0]}_{srt_acq.split('.')[1]}_{sat}_FRP.png")

def ACQ_datetime(date_want, sat):
    """ Calculates acquistation datetime for filename """
    
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
    """Calculates the day of year from user input year, month,day. """

    doy = date(int(yr), int(mn), int(dy)).timetuple().tm_yday
    return doy

#%% User Defined Variables
if __name__ == "__main__":
    t_start = datetime.now()
    save_dir = 'C:/Users/k2262276/Documents/FENIX_SAT/'
    plot_dir = 'C:/Users/k2262276/Documents/FENIX_Plot/'
    
    day_range = 7
    daynight = 'D' # 'D'/'N'
    satellites = ['NOAA-20', 'SNPP','AQUA','TERRA','LANDSAT'] # 'AQUA' / 'TERRA' / 'SNPP' / 'NOAA-20' / 'LANDSAT'
    pickle_lake = [-90.18333, 51.466667] # Longitude, Latitude
    date_want = datetime(2023,6,28) # datetime.now() 
    
    #%% FIRMS Scraping and Plotting
    for sat in satellites:
        df = data_retrieve(sat, pickle_lake, daynight, day_range, date_want)   
        rel_orb, un = unique_orbit(df)
        
        for ro_val in np.arange(np.max(rel_orb)+1):
            if len(un) > 0:
                srt = np.datetime64(un[[i for i, j in enumerate(rel_orb) if j == ro_val]].min(),'m')
                end = np.datetime64(un[[i for i, j in enumerate(rel_orb) if j == ro_val]].max(),'m')
                rel_orb_df = df[df['acq_datetime'].isin(un[[i for i, j in enumerate(rel_orb) if j == ro_val]])]        
                plotting(rel_orb_df, pickle_lake, srt, end, sat, plot_dir)
            else:
                print(f"{sat} detects no fires within 250km of {pickle_lake[1]}, {pickle_lake[0]}, in the {day_range} days leading up to {date_want.strftime('%D')}.")
    