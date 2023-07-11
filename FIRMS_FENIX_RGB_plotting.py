# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:49:07 2023

@author: k2262276
"""
from cartopy.geodesic import Geodesic
from bs4 import BeautifulSoup # Library to scrap information from online
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from datetime import date, datetime, timedelta
from tqdm import tqdm # Library to create progress bars for loops/functions
from glob import glob
from skimage import exposure
from pyhdf.SD import SD, SDC

import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.ticker as mticker
import numpy as np
import netCDF4 as nc

import warnings, requests, subprocess, os, calendar, sys

#%% Preamble
warnings.filterwarnings('ignore')
plt.close('all')
startTime = datetime.now()

#%% UDF
def RGB_Bands(sat):
    if sat.upper() == 'SNPP' or sat.upper() == 'NOAA-20':
        Red = 'M05'
        Green = 'M04'
        Blue = 'M03'
    else:
        Red = '1'
        Green = '4'
        Blue = '3'
    return Red, Green, Blue

def date_finder(yr, doy):
    """
    Calculates the month and day from user input of year and day of year.

    Parameters
    ----------
    yr : Float
        Year number.
    doy : interger
        Day of Year.

    Returns
    -------
    M_NUM : String
        Month Number.
    DAY : String
        Day Number.

    """
    year = float(yr)
    doy = int(doy)
    
    # Is the year a leap year?
    if calendar.isleap(year):
        Y = 1
        
    else:
        Y = 0
    
    # Check if doy is valid
    if Y == 0 and (doy > 365 or doy < 1):
        sys.exit(str(int(year)) + " is not a leap year. Doy value is invalid (> 365 / < 1)")
    elif Y == 1 and (doy > 366 or doy < 1):
        sys.exit(str(int(year)) + " is a leap year. Doy value is invalid (> 366 / < 1)")

    
    Mon_Len = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # mon len=  1   2   3   4   5   6   7   8   9  10  11  12

    if Y == 1:
        Mon_Len[1] += 1
        
        # determine the number of the day in the month
    for N, M in enumerate(Mon_Len):
        if doy > M:
            doy = doy-M
        else:
            M_NUM = N+1
            break

    if len(str(doy)) == 1:
        DAY = '0' + str(doy)
    else:
        DAY = str(doy)

    if len(str(M_NUM)) == 1:
        M_NUM = '0' + str(M_NUM)
    else:
        M_NUM = str(M_NUM)

    return M_NUM, DAY

def glob_checker(directory, search_string):
    """
    Creates list of files of certain file extension from directory.

    Parameters
    ----------
    directory : String
        Absolute Directory to be checked.
    search_string : String
        String to use with wildcards to search within the directory.

    Returns
    -------
    abs_fname_list : List of Strings
        List of strings of files within the directory that meets the search 
        criteria.
    """
    to_check = directory + search_string
    rel_fname = list(f.split('\\')[-1] for f in glob(to_check))
    abs_fname_list = []
    for fname in rel_fname:
        abs_fname_list.append(os.path.join(directory,fname))
    return abs_fname_list

def dist_det(north, east, west, south):
    dist = np.ceil(north)-north
    if np.ceil(west)-west > dist:
        dist = np.ceil(west)-west
    if south - np.floor(south) > dist:
        dist = south - np.floor(south)
    if east - np.floor(east) > dist:
        dist = east - np.floor(east)
    
    return dist

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

def circle_calc(pickle_lake, radius=200):
    gd = Geodesic()
    cp = gd.circle(lon=pickle_lake[0], lat=pickle_lake[1], radius = radius*1e3)
    cp = np.vstack((cp, cp[0,:]))
    return cp

def plotting_FRP(rel_orb_df, north, east, west, south, pickle_lake, red_lake, radius, srt, end, sat, plot_dir):
    #dist = dist_det(north, east, west, south)
    cp = circle_calc(pickle_lake, radius)

    fig, ax = plt.subplots(figsize=(15,10), subplot_kw={"projection": ccrs.Mercator()})
    ax.set_extent([np.floor(west*2)/2, np.ceil(east*2)/2, np.floor(south*2)/2, np.ceil(north*2)/2], crs=ccrs.PlateCarree())
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
    plt.plot(cp[:,0], cp[:,1],
             transform=ccrs.PlateCarree(),
             zorder = 0)
    cbar = plt.colorbar(sca)
    cbar.set_label(label = '$log_{10}(frp)$', 
                  size = 14, weight = 'bold')

    ax.add_feature(cf.BORDERS, linestyle=':')
    ax.add_feature(cf.COASTLINE)
    ax.add_feature(cf.STATES, linestyle='-.')
    ax.add_feature(cf.LAKES, alpha=0.5)
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
    gl.xlocator = mticker.FixedLocator(np.arange(np.floor(west*2)/2,
                                                 np.ceil(east*2)/2,
                                                 ((np.ceil(east*2)/2) - (np.floor(west*2)/2))/5))
    gl.ylocator = mticker.FixedLocator(np.arange(np.floor(south*2)/2,
                                                 np.ceil(north*2)/2,
                                                 ((np.ceil(north*2)/2) - (np.floor(south*2)/2))/5))
    
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    lon_formatter = LongitudeFormatter(number_format='.2f',
                                   degree_symbol='',
                                   dateline_direction_label=True)
    lat_formatter = LatitudeFormatter(number_format='.2f',
                                  degree_symbol='')
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    plt.title(f"{sat} | Radius: {radius}km\n{srt} --> {end} UTC")
    
    srt_acq = ACQ_datetime(srt, sat)
    end_acq = ACQ_datetime(end, sat)
    
    plt.savefig(f"{plot_dir}{sat}_{srt_acq.split('.')[0]}_{srt_acq.split('.')[1]}_{end_acq.split('.')[1]}_FRP.png", bbox_inches='tight', dpi=300, pad_inches = 0)
    plt.close()
    
def doy_finder(yr, mn, dy):
    """Calculates the day of year from user input year, month,day.
    Input: yr: year (yyyy)
    Input: mn: month (mm)
    Input: dy: day (dd)

    Output: doy: Day of Year (1st Jan = 1)"""

    doy = date(int(yr), int(mn), int(dy)).timetuple().tm_yday
    return doy

def ACQ_datetime(date_want, sat):
    year_wnt = str(date_want).split('-')[0]
    month_wnt = str(date_want).split('-')[1]
    day_wnt = str(date_want).split('-')[2].split('T')[0]
    hour_wnt = str(date_want).split('T')[1].split(':')[0]
    minute_wnt = str(date_want).split('T')[1].split(':')[1]

    doy = doy_finder(year_wnt, month_wnt, day_wnt)

    if sat.upper() == 'SNPP' or sat.upper() == 'NOAA-20':
        min_floor = str((int(minute_wnt) // 6)*6)
    else:
        min_floor = str((int(minute_wnt) // 5)*5)

    if len(min_floor) != 2:
        min_floor = '0' + min_floor
        
    time_wnt = hour_wnt + min_floor
    
    acq_datetime = 'A' + str(year_wnt) + str(doy) + "." + time_wnt
    
    return acq_datetime

def listFD(url, ext=''):
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    return [url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]

def runcmd(cmd, verbose = False, *args, **kwargs):

    process = subprocess.Popen(
        cmd,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        text = True,
        shell = True
    )
    std_out, std_err = process.communicate()
    if verbose:
        print(std_out.strip(), std_err)
    pass

def Product_download(collection, product, year, doy, start_time, end_time, ext, save_dir, token):    
    
    master_url = 'https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/allData/'
    search_url = master_url + collection + '/' + product + '/' + year + '/' + doy + '/'

    files = []
    for file in listFD(search_url, ext):
        #for swath in file:
        if product[0].upper() == 'V':
            acq_time = file.split('.')[-4]
        else:
            acq_time = file.split('.')[-5]
        if (int(acq_time) >= int(start_time)) and (int(acq_time) <= int(end_time)):
               files.append(file)
    
    files = list(set(files))
    for file in tqdm(files, unit='files', bar_format="{desc}Downloading:{percentage:3.0f}%|{bar}|{n_fmt}/{total_fmt} [{elapsed}<{remaining}]"):
    #for file in files:
        url = file.split('///')[0] + '/' + file.split('///')[-1]
         
        if not os.path.isfile(save_dir + url.split('/')[-1]):
            runcmd(f'wget -e robots=off -m -np -nv -nd -R .html,.tmp -nH --cut-dirs=3 {url} --header "Authorization: Bearer {token}" -P {save_dir}', verbose = True)    

def tens_str(num):
    if num >= 10:
        num_str = str(num)
    else:
        num_str = '0'+str(num)
    return num_str

def data_download(file_id, sat, save_dir):
    if sat.upper() == 'AQUA' or sat.upper() == 'TERRA':
        collection = '61'
        ext = '.hdf'
        if sat.upper() == 'TERRA':
            prefix = 'MOD'
        else:
            prefix = 'MYD'
    else:
        collection = '5200'
        ext = '.nc'
        if sat.upper() == 'SNPP':
            prefix = 'VNP'
        else:
            prefix = 'VJ1'
    
    m_num, d_num = date_finder(file_id[1:5], file_id[5:8])
    srt_datetime = datetime(int(file_id[1:5]), int(m_num), int(d_num), int(file_id[9:11]), int(file_id[11:13])) - timedelta(minutes=1)
    end_datetime = datetime(int(file_id[1:5]), int(m_num), int(d_num), int(file_id[9:11]), int(file_id[11:13])) + timedelta(minutes=1)
    
    srt_time = tens_str(srt_datetime.hour) + tens_str(srt_datetime.minute)
    end_time = tens_str(end_datetime.hour) + tens_str(end_datetime.minute)

    if sat.upper() == 'SNPP' or sat.upper() == 'NOAA-20':
        Product_download(collection, prefix+'03MOD_NRT', file_id[1:5], file_id[5:8], srt_time, end_time, ext, save_dir, token)
        Product_download(collection, prefix+'02MOD_NRT', file_id[1:5], file_id[5:8], srt_time, end_time, ext, save_dir, token)
        Product_download(collection, prefix+'03IMG_NRT', file_id[1:5], file_id[5:8], srt_time, end_time, ext, save_dir, token)
        Product_download(collection, prefix+'02IMG_NRT', file_id[1:5], file_id[5:8], srt_time, end_time, ext, save_dir, token)  
    elif sat.upper() == 'AQUA' or sat.upper() == 'TERRA':
        Product_download(collection, prefix+'021KM', file_id[1:5], file_id[5:8], srt_time, end_time, ext, save_dir, token)
        Product_download(collection, prefix+'03', file_id[1:5], file_id[5:8], srt_time, end_time, ext, save_dir, token)

def kron(array,n):
    """
    expanding the 2D array, to match differing resolutions

    Parameters
    ----------
    array : Array
        Array to be expanded.
    n : Integar
        expansion factor (2 = double both sides).

    Returns
    -------
    data : Array
        Expanded Array.

    """
    
    data = np.kron(array, np.ones((n,n)))
    return data

def VIIRS_M_Band(mref_file, scle, Red='M05', Green='M04', Blue='M03'):
    """
    Parameters
    ----------
    mref_file : String
        Moderate-resolution Calibrated Radiances.
    scle : String
        Toggle for scaling array ('scaled' allows for scaling via percentage, 
                                  otherwise it is MinMax).

    Returns
    -------
    M03 : Array
        Array of the MO3 Band (Blue), at Image-Resolution.
    M04 : Array
        Array of the M04 Band (Green), at Image-Resolution.
    M05 : Array
        Array of the M05 Band (Red), at Image-Resolution.

    """
    if scle.lower() == 'scaled':
        B = kron(scalePercent(VIIRS_Band_nc(mref_file, Blue),2),2)
        G = kron(scalePercent(VIIRS_Band_nc(mref_file, Green),2),2)
        R = kron(scalePercent(VIIRS_Band_nc(mref_file, Red),2),2)
    else:
        B = kron(scaleMinMax(VIIRS_Band_nc(mref_file, Blue)),2)
        G = kron(scaleMinMax(VIIRS_Band_nc(mref_file, Green)),2)
        R = kron(scaleMinMax(VIIRS_Band_nc(mref_file, Red)),2)

    return B, G, R

def VIIRS_Band_nc(abs_filename, band_name):
    """
    Loading in VIIRS Spectral Radiance Data for both VNP02MOD/IMG 

    Parameters
    ----------
    abs_filename : String
        Absolute Filename for the VIIRS Calibrated Radiance files.
    band_name : String
        Band Name to import.

    Returns
    -------
    scaled_band : Array
        Scaled and Offset Calibrated Radiance data.

    """
    ds = nc.Dataset(abs_filename)
    Raw_Band = ds['observation_data'][band_name][:]
    offset = ds['observation_data'][band_name].add_offset
    scale_factor = ds['observation_data'][band_name].scale_factor
    scaled_band = (Raw_Band/scale_factor) + offset
    flags = ds['observation_data'][band_name].flag_values
    contains = np.in1d(scaled_band,flags).reshape((scaled_band.shape[0],scaled_band.shape[1]))
    scaled_band[contains == 1] = np.nan
    
    ds.close()
    return scaled_band

def scaleMinMax(x):
    """
    Produces a Scaling Index between the smallest and largest values


    Parameters
    ----------
    x : Array
        Array to be scaled.

    Returns
    -------
    scle : Array
        Scaled Array.

    """
    scle = (x - np.nanmin(x))/(np.nanmax(x) - np.nanmin(x))
    return scle

def scalePercent(x, pcnt):
    """
    Produces a Scaling Index between the pcnt% and (100-pcnt)% values

    Parameters
    ----------
    x : Array to be scaled
        DESCRIPTION.
    pcnt : Integar
        How much to trim from the array.

    Returns
    -------
    scle : Array
        Scaled .

    """
    scle = (x - np.nanpercentile(x, pcnt)) / (np.nanpercentile(x, 100 - pcnt) - np.nanpercentile(x, pcnt))
    
    scle[scle >= 1] = 1
    scle[scle <= 0] = 0
    
    return scle

def VIIRS_I_Band(igeo_file, iref_file):
    """
    Loads the data from the Image-Resolution Files

    Parameters
    ----------
    igeo_file : String
        Absolute Filepath of the Image Geolocation data.
    iref_file : String
        Absolute Filepath of the Image Radiances data.

    Returns
    -------
    ilat : Array
        Image-resolution Latitude data.
    ilon : Array
        Image-resolution Longitude data.
    I01 : Array
        Array of I01 Data (Red).
    I02 : Array
        Array of I02 Data (NIR).

    """
    
    ilat, ilon = GEO_nc(igeo_file)
    I01 = scaleMinMax(VIIRS_Band_nc(iref_file, 'I01'))
    I02 = scaleMinMax(VIIRS_Band_nc(iref_file, 'I02'))
    return ilat, ilon, I01, I02

def GEO_nc(abs_filename):
    """
    Loading in VIIRS Geolocation Data for both VNP03MOD/IMG 
    lat, lon = GEO_nc(abs_filename)

    Parameters
    ----------
    abs_filename : String
        Absolute Filename of Geolocation data.

    Returns
    -------
    lat : Array
        Array of Latitudes.
    lon : Array
        Array of Longitudes.

    """
    ds = nc.Dataset(abs_filename)
    lat = ds['geolocation_data']['latitude'][:]
    lon = ds['geolocation_data']['longitude'][:]
    lon[lon < 0] += 360 # Re-aligning longitudes due to +180/-180 issue
    ds.close()
    return lat, lon

def Colour_Composite(files, Red, Green, Blue, sat):
    if sat.upper() == 'SNPP' or sat.upper() == 'NOAA-20':
        B, G, R = VIIRS_M_Band(files[1], 'full', Red, Green, Blue)
        ilat, ilon, __, __ = VIIRS_I_Band(files[2], files[0])
        
    else:
        R = MODIS_Calbrated_Band(files[0], int(Red))
        G = MODIS_Calbrated_Band(files[0], int(Green))
        B = MODIS_Calbrated_Band(files[0], int(Blue))
        ilat, ilon = MODIS_Geolocation(files[1])

    RGB = np.dstack((R.data, G.data, B.data))
    p2, p98 = np.nanpercentile(RGB, (2, 98))
    rgbStretched = exposure.rescale_intensity(RGB, in_range=(p2, p98))
    rgbStretched = exposure.adjust_gamma(rgbStretched, 0.5)
    #RGB_mask = np.dstack((R.mask, G.mask, B.mask))
    #RGBA = np.dstack([np.floor(rgbStretched*255), ((1.-np.any(RGB_mask, axis=2))*255).astype(int)])
    colorTuple = rgbStretched.reshape((rgbStretched.shape[0]*rgbStretched.shape[1],rgbStretched.shape[2]))
    print(f"Making RGB VIIRS DATA: {datetime.now() - startTime}")
    return rgbStretched, colorTuple, ilat, ilon

def file_downloader(un_time, un, sat, clip_df, sat_pre):
    file_id = ACQ_datetime(un[un_time], sat)
    rel_orb_df = clip_df[clip_df['acq_datetime'].isin(un[[i for i, j in enumerate(rel_orb) if j == rel_orb[un_time]]])]        

    files = glob_checker(save_dir, f'{sat_pre}*{file_id}*')
    if not files:
        data_download(file_id, sat, save_dir)
        print(f"Downloading VIIRS DATA: {datetime.now() - startTime}")

        files = glob_checker(save_dir, f'*{file_id}*')
    
    return files, rel_orb_df, file_id

def plotting_FCC(north, east, west, south, ilon, ilat, rgbStretched, colorTuple, rel_orb_df, pickle_lake, red_lake,file_id, Red, Green, Blue, plot_dir, sat_pre):

    dist = dist_det(north, east, west, south)

    fig, ax = plt.subplots(figsize=(15,10), subplot_kw={"projection": ccrs.Mercator()})
    ax.set_extent([west-dist, east+dist, south-dist, north+dist], crs=ccrs.PlateCarree())
    plt.pcolormesh(ilon, ilat, rgbStretched[:,:,0], color = colorTuple, transform=ccrs.PlateCarree())#, linewidth=0)
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
    ax.add_feature(cf.BORDERS, edgecolor='cyan')
    ax.add_feature(cf.COASTLINE, edgecolor='cyan')
    ax.add_feature(cf.STATES, edgecolor='cyan')
    ax.set_facecolor("gainsboro")
    ax.set_title(f'\n{file_id} - {Red}/{Green}/{Blue}')

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
    gl.xlocator = mticker.FixedLocator(np.arange(np.floor(west*2)/2-dist,
                                                 np.ceil(east*2)/2+(dist*2),
                                                 ((np.ceil(east*2)/2+dist) - (np.floor(west*2)/2-(dist*2)))/5))
    gl.ylocator = mticker.FixedLocator(np.arange(np.floor(south*2)/2-dist,
                                                 np.ceil(north*2)/2+(dist*2),
                                                 ((np.ceil(north*2)/2+dist) - (np.floor(south*2)/2-(dist*2)))/5))
    
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    lon_formatter = LongitudeFormatter(number_format='.2f',
                                   degree_symbol='',
                                   dateline_direction_label=True)
    lat_formatter = LatitudeFormatter(number_format='.2f',
                                  degree_symbol='')
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    plt.savefig(f"{plot_dir}{sat}_{file_id.split('.')[0]}_{file_id.split('.')[1]}_{Red}{Green}{Blue}.png", bbox_inches='tight', dpi=300, pad_inches = 0)
    plt.close()

def MODIS_Band_Calibration(data, ch_index):
    scales=data.attributes()['radiance_scales']
    offsets=data.attributes()['radiance_offsets']
    ch_scale=scales[ch_index]
    ch_offset=offsets[ch_index]
    ch_calibrated =(data[ch_index,:,:] - ch_offset)*ch_scale
    return ch_calibrated

def MODIS_Calbrated_Band(FILE_NAME, BAND_NUM):
    hdf = SD(FILE_NAME, SDC.READ)
    
    if BAND_NUM <= 2:
        data_var_name = 'EV_250_Aggr1km_RefSB'
        bands_var_name = hdf.select('Band_250M')
    elif BAND_NUM <= 7:
        data_var_name = 'EV_500_Aggr1km_RefSB'
        bands_var_name = hdf.select('Band_500M')
    else:
        data_var_name = 'EV_1KM_RefSB'
        bands_var_name = hdf.select('Band_1KM_RefSB')
        
    band_nums = bands_var_name.get()
    ch_index=np.searchsorted(band_nums,BAND_NUM)
    ch_index = int(ch_index)

    band_data = hdf.select(data_var_name) # select sds
    ch_calibrated = MODIS_Band_Calibration(band_data, ch_index)
    hdf.end()
    return ch_calibrated

def MODIS_Geolocation(GEO_FILE_NAME):
    hdf_geo = SD(GEO_FILE_NAME, SDC.READ)
    
    # Read geolocation dataset from MOD03 product.
    lat = hdf_geo.select('Latitude')
    latitude = lat[:,:]
    lon = hdf_geo.select('Longitude')
    longitude = lon[:,:]
    hdf_geo.end()
    
    return latitude, longitude

#%% User Derived Varaibles
token = 'eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJBUFMgT0F1dGgyIEF1dGhlbnRpY2F0b3IiLCJpYXQiOjE2ODcyNjQ4NzMsIm5iZiI6MTY4NzI2NDg3MywiZXhwIjoxNzAyODE2ODczLCJ1aWQiOiJ3bWFzbGFua2EiLCJlbWFpbF9hZGRyZXNzIjoid2lsbGlhbS5tYXNsYW5rYUBrY2wuYWMudWsiLCJ0b2tlbkNyZWF0b3IiOiJ3bWFzbGFua2EifQ.SpghatLJR5t7s14Fw93sVT-HGLcLWceQvDS1HPc93Gc'
save_dir = 'C:/Users/k2262276/Documents/FENIX_SAT/'
plot_dir = 'C:/Users/k2262276/Documents/FENIX_Plot/'

north = -15.6097033
south = -22.4241096

west = 25.2373
east = 33.0683413
daynight = 'D' # 'D'/'N'
day_range = 2
pickle_lake = [-90.18333, 51.466667]
red_lake = [-93.793056, 51.0679222]
radius = 0.1 # km 

sat_all = ['NOAA-20']#['TERRA', 'AQUA', 'SNPP', 'NOAA-20']

for sat in sat_all:
    Red, Green, Blue = RGB_Bands(sat)
    
    if sat.upper() == 'SNPP':
        sat_pre = 'VNP'
    elif sat.upper() == 'NOAA-20':
        sat_pre = 'VJ1'
    elif sat.upper() == 'TERRA':
        sat_pre = 'MOD'
    elif sat.upper() == 'AQUA':
        sat_pre = 'MYD'
    
    #%% Retrieved FRP data from FIRMS
    df = data_retrieve(sat, north, east, west, south, daynight, day_range)   
    clip_df = clipping_df(df)    
    rel_orb, un = unique_orbit(clip_df)
    print(f"Calculating Relative Orbit | {sat}: {datetime.now() - startTime}")
    
    #%% Plotting the FRP data from FIRMS
    for ro_val in np.arange(np.max(rel_orb)+1):
        srt = np.datetime64(un[[i for i, j in enumerate(rel_orb) if j == ro_val]].min(),'m')
        end = np.datetime64(un[[i for i, j in enumerate(rel_orb) if j == ro_val]].max(),'m')
    
        rel_orb_df = clip_df[clip_df['acq_datetime'].isin(un[[i for i, j in enumerate(rel_orb) if j == ro_val]])]        
        plotting_FRP(rel_orb_df, north, east, west, south, pickle_lake, red_lake, radius, srt, end, sat, plot_dir) 
    
    
    #%% Downloading and Plotting the Colour Composite
    for un_time in np.arange(len(un)):
        files, rel_orb_df, file_id = file_downloader(un_time, un, sat, clip_df, sat_pre)
        rgbStretched, colorTuple, ilat, ilon = Colour_Composite(files, Red, Green, Blue, sat)
        plotting_FCC(north, east, west, south, ilon, ilat, rgbStretched, colorTuple, rel_orb_df, pickle_lake, red_lake,file_id, Red, Green, Blue, plot_dir, sat)
        print(f"Saved RGB Image | {sat}: {datetime.now() - startTime}")
    