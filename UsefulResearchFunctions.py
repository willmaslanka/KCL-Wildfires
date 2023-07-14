# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 17:34:01 2023
@author: Will Maslanka

Useful Research Functions
"""

# Need to get water vapour units
# Need to finish comments

#%% Modules
import os, subprocess, requests, sys, calendar
import s3fs # Module to interface with Amazon Simple Storage Service (S3)
from glob import glob
from sentinelsat import SentinelAPI # Library to access Sentinel API
from bs4 import BeautifulSoup # Library to scrap information from online
from tqdm import tqdm # Library to create progress bars for loops/functions
from packaging.version import parse # Library to access core utilities for Python packages
import numpy as np
import datetime
from pyproj import Proj

#%% User Defined Functions
########## Misc. Functions ##########
def glob_checker(directory, search_string):
    """Creates list of files of certain file extension from directory"""
    to_check = directory + search_string
    #glob_check = glob(to_check)
    rel_fname = list(f.split('\\')[-1] for f in glob(to_check))
    abs_fname_list = []
    for fname in rel_fname:
        abs_fname_list.append(os.path.join(directory,fname))
    return abs_fname_list

def closest_value(input_list, input_value):
    """ Returns the index of the closest value to a user submitted value """
    arr = np.asarray(input_list)
    i = (np.abs(arr - input_value)).argmin()
    return i

def tens_str(num):
    """ Returns a string with leading zeros.
    Works for positive numbers and zero. """
    if num >= 10:
        num_str = str(num)
    else:
        num_str = '0'+str(num)
    return num_str

######## Sentinel-5P Tropomi Download ########
def download_Sentinel5P(bbox, date_trop_start, date_trop_end, start_time_T, end_time_T, S5P_dir):
    """
    Master Function to download Sentinel-5P data to local drive

    Parameters
    ----------
    bbox : list of string
        Bounding Box of interest, [west, south, east, north].
    date_trop_start : String
         Date before Acquistion in Sentinel-5P naming convention.
    date_trop_end : String
         Date of Acquistion in Sentinel-5P naming convention.
    start_time_T : String
        Start time of Tropomi Acqusition.
    end_time_T : String
        End time of Tropomi Acqusition.
    S5P_dir : string
        Directory to save the data in.

    Returns
    -------
    None.

    """
    product_abbreviation = get_tropomi_product_abbreviation('CO')
    get_tropomi_files(bbox[0], bbox[2], bbox[1], bbox[3],
                      date_trop_start, date_trop_end, start_time_T, end_time_T,
                      product_abbreviation, 'Reprocessing', S5P_dir)
    print('Sentinel-5P Downloaded')
    
def get_tropomi_product_abbreviation(product):
    """Gets Tropomi Products"""
    if product == 'CO':
        product_abbreviation = 'L2__CO____'
    elif product == 'NO2':
        product_abbreviation = 'L2__NO2___'
    elif product == 'SO2':
        product_abbreviation = 'L2__SO2___'
    elif product == 'HCHO':
        product_abbreviation = 'L2__HCHO__'
    elif product == 'AI':
        product_abbreviation = 'L2__AER_AI'
    elif product == 'ALH':
        product_abbreviation = 'L2__AER_LH'
    return product_abbreviation

def get_tropomi_files(west_lon, east_lon, south_lat, north_lat, start_date, end_date, start_time, end_time, product_abbreviation, latency, save_path):
    """
    Function to Download Tropomi Files

    Parameters
    ----------
    west_lon : string
        West-most longitude.
    east_lon : string
        East-most longitude.
    south_lat : string
        Southern-most latitude.
    north_lat : string
        Northern-most string
    date_trop_start : String
         Date before Acquistion in Sentinel-5P naming convention.
    date_trop_end : String
         Date of Acquistion in Sentinel-5P naming convention.
    start_time_T : String
        Start time of Tropomi Acqusition.
    end_time_T : String
        End time of Tropomi Acqusition.
    product_abbreviation : string
        Trace Gas Wanted. Look in get_tropomi_product_abbreviation
    latency : string
        Product of Sentinel-5P - 'Reprocessing' get reprocessed data.
    save_path : string
        Directory to save the data in.

    Returns
    -------
    None.

    """
    # Query S5P Data Hub and list file names matching user-entered info
    file_name_list, file_size_list, products = tropomi_list_files(west_lon, east_lon, south_lat, north_lat, start_date, end_date, start_time, end_time, product_abbreviation, latency)
    
    # Print list of available file names/sizes
    if len(file_name_list) > 0:
        print('\nList of available data files (file size):')
        for file, size in zip(file_name_list, file_size_list):
            print(file, ' (', size, ')', sep='')
        
        # Print directory where files will be saved, then download
        print('\nData files will be saved to:', save_path)
        tropomi_download_files(products, save_path)
    else:
        print('\nNo files retrieved.  Check settings and try again.')
        
def tropomi_list_files(west_lon, east_lon, south_lat, north_lat, start_date, end_date, start_time, end_time, product_abbreviation, latency):
    """
    Creates list of data of Sentinel-5P files to download

    Parameters
    ----------
    west_lon : string
        West-most longitude.
    east_lon : string
        East-most longitude.
    south_lat : string
        Southern-most latitude.
    north_lat : string
        Northern-most string
    date_trop_start : String
         Date before Acquistion in Sentinel-5P naming convention.
    date_trop_end : String
         Date of Acquistion in Sentinel-5P naming convention.
    start_time_T : String
        Start time of Tropomi Acqusition.
    end_time_T : String
        End time of Tropomi Acqusition.
    product_abbreviation : string
        Trace Gas Wanted. Look in get_tropomi_product_abbreviation
    latency : string
        Product of Sentinel-5P - 'Reprocessing' get reprocessed data.

    Returns
    -------
    file_name_list : list of string
        List of filenames.
    file_size_list : list of string
        List of filesizes.
    products : TYPE
        DESCRIPTION.

    """
    # Access S5P Data Hub using guest login credentials
    api = SentinelAPI('s5pguest', 's5pguest', 'https://s5phub.copernicus.eu/dhus')
   
    # Query API for specified region, start/end dates, data product
    footprint = 'POLYGON((' + west_lon + ' ' + south_lat + ',' + east_lon + ' ' + south_lat + ',' + east_lon + ' ' + north_lat + ',' + west_lon + ' ' + north_lat + ',' + west_lon + ' ' + south_lat + '))'
    try:
        products = api.query(area=footprint, date=(start_date + 'T' + start_time, end_date + 'T' + end_time), producttype=product_abbreviation, processingmode=latency)
    except:
        print('Error connecting to SciHub server. This happens periodically. Run code again.')
    
    # Convert query output to pandas dataframe (df) (part of Sentinelsat library)
    products_df = api.to_dataframe(products)
    
    # Extract data file names from dataframe to list
    if len(products_df) > 0:
        file_name_list = products_df['filename'].tolist()
        file_size_list = products_df['size'].tolist()
    else:
        file_name_list = []
        file_size_list = []
    
    return file_name_list, file_size_list, products

def tropomi_download_files(products, save_path):
    """Downloads Tropomi Files"""
    # Query S5P Data Hub using guest login credentials
    api = SentinelAPI('s5pguest', 's5pguest', 'https://s5phub.copernicus.eu/dhus')
    
    # Download data files to specified subdirectory
    # Note: Sentinelsat library includes tqdm download progress bar
    try:
        api.download_all(products, save_path)
    except KeyboardInterrupt:
        print('\nDownload was interrupted by user.')
        
######## LAADS Download ######## 
def Product_download(collection, product, year, doy, start_time, end_time, ext, save_dir, token):
    """JPSS Function to download data"""
    master_url = 'http://ladsweb.modaps.eosdis.nasa.gov/archive/allData/'
    search_url = master_url + collection + '/' + product + '/' + year + '/' + doy + '/'

    files = []
    for file in listFD(search_url, ext):
        #for swath in file:
        acq_time = file.split('///')[-1].split('.')[2]
        if (int(acq_time) >= int(start_time)) and (int(acq_time) <= int(end_time)):
               files.append(file)
    
    files = list(set(files))
    for file in tqdm(files, unit='files', bar_format="{desc}Downloading:{percentage:3.0f}%|{bar}|{n_fmt}/{total_fmt} [{elapsed}<{remaining}]"):
    #for file in files:
        url = file.split('///')[0] + '/' + file.split('///')[-1]
         
        if not os.path.isfile(save_dir + url.split('/')[-1]):
            runcmd(f'wget -e robots=off -m -np -q -nd -R .html,.tmp -nH --cut-dirs=3 {url} --header "Authorization: Bearer {token}" -P {save_dir}', verbose = True)    

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

def NRT_Product_download(collection, product, year, doy, start_time, end_time, ext, save_dir, token):    
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

########### AWS VIIRS AOD EPS ##########
def get_product_path_AODEPS(product, processing, satellite):
    """ Function to create product string """
    # Set AWS NODD product name
    if product == 'AOD':
        product_name = 'Aerosol_Optical_Depth'
    elif product == 'ADP':
        product_name = 'Aerosol_Detection'
        
    # Set AWS NODD processing suffix
    if processing == 'reprocessed':
        suffix = '_EDR_Reprocessed'
    elif processing == 'operational':
        suffix = '_EDR'
    
    # Create AWS NOD product path name
    product_path = satellite + '_VIIRS_' + product_name + suffix

    return product_path

def aws_viirs_list(year, month, day, start_hour, start_min, end_hour, end_min, satellite, product, processing):
    """ Function to get list of VIIRS data from AWS server for the given time period.
    yyyy/mm/dd start time (hh:mm) to end time (hh:mm)
    for a given satellite, product, and processing level."""
    
    # Construct aerosol product path for AWS NODD
    product_path = get_product_path_AODEPS(product, processing, satellite)

    # Access AWS using anonymous credentials
    aws = s3fs.S3FileSystem(anon=True)
    
    # Create list of file names for entire day (~550)
    day_files = aws.ls('noaa-jpss/' + satellite + '/VIIRS/' + product_path + '/' + year + '/' + month + '/' + day + '/', refresh=True)  
    
    # Create list of subsetted file names that fall within specified time period(s)
    data = []
    for file in day_files:
        file_time = file.split('_')[-3][9:13]
        if file_time >= (start_hour + start_min) and file_time <= (end_hour + end_min):
            data.append(file)
        
    return data

def get_viirs_files(year, month, day, start_hour, start_min, end_hour, end_min, satellite, product, processing, save_path):
    """ Downloads the data from AWS"""
    # Query AWS VIIRS archive and print names/sizes of available aerosol files
    data = aws_viirs_list(year, month, day, start_hour, start_min, end_hour, end_min, satellite, product, processing)
    al_down = glob_checker(save_path,f'*{product}*.nc')
    
    data_download = []
    for en_a, file in enumerate(data):
        splitter = file.split('/')[-1]
        abs_split = save_path + splitter
        if abs_split not in al_down:
            data_download.append(file)

    if len(data_download) > 0:   
        # Access AWS using anonymous credentials
        aws = s3fs.S3FileSystem(anon=True)
        tot_file_size = 0;
        # Print list of available data files
        for file in data_download:
            file_size = aws.size(file)
            # sep='' removes extra spaces b/w print elements
            #print(file.split('/')[-1], ' (', int(file_size/1.0E6), ' MB)', sep='')
            tot_file_size = tot_file_size + int(file_size/1.0E6)
        
        # Print directory where files will be saved
        print('Data files will be saved to: ' + str(save_path))
        
        # Ask user if they want to download the available data files
        # If yes, download files to specified directory
        download_question = f'Downloading the {str(len(data_download))} files? ({str(np.round(tot_file_size/1000,2))} GB) for {day}/{month}/{year}.'# Would you like to download the ' + str(len(data)) + ' files? (' + str(np.round(tot_file_size/1000,2)) + ' GB)\nType "yes" or "no" and hit "Enter"\n'
        print(download_question)
        download_files = 'yes' #input(download_question)
        if download_files in ['yes', 'YES', 'Yes', 'y', 'Y']:
            
            # Display progress bar using tqdm library
            # Flush buffer if Python version < v3.9 to avoid glitch in tqdm library
            if parse(sys.version.split(' ')[0]) < parse('3.9'):
                sys.stdout.flush()
            for name in tqdm(data_download, unit='files', bar_format="{desc}Downloading:{percentage:3.0f}%|{bar}|{n_fmt}/{total_fmt} [{elapsed}<{remaining}]"):
                # Set save_path + file_name as pathlib.Path object and convert to string
                #full_path = str(save_path / name.split('/')[-1])
                full_path = save_path + str(name.split('/')[-1])
                # Download file from AWS archive
                aws.get(name, full_path)
            print('Download complete!')
        else:
            print('Files are not being downloaded.')
    else:
        print('No files retrieved. Check settings and try again.')

########### FRP Atmospheric Correction ##########
def get_Atmostrans_VIIRS(vangle,water_vapour):
    """
    Look Up Table to calculate the atmospheric transmissivity for VIIRS, using
    an user inputted angle, and water vapour.

    Parameters
    ----------
    vangle : float
        Satellite View Angle.
    water_vapour : float
        Water vapour in kg /m2 units.

    Returns
    -------
    atmos_trans : float
        Atmospheric Transmissivity.
    """
    H2O_list=[5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
    TAU_list=[0.118223, 0.143704, 0.166257, 0.186838, 0.206053, 0.224324, 0.244870, 0.266351, 0.288009, 0.305552, 0.324089, 0.339992]
    P0_list =[+0.015990081, +0.017316598, +0.018890205, +0.020397422, +0.021573963, +0.022437754, +0.023041454, +0.023512723, +0.023815163, +0.024093296, +0.024290476, +0.024453227]
    P1_list =[+1.125117200, +1.104931300, +1.081749500, +1.060689700, +1.042883900, +1.028166400, +1.017686100, +1.009994700, +1.004233400, +0.997618210, +0.992723030, +0.987414580]
    P2_list =[-0.080837809, -0.073449551, -0.063037274, -0.052937644, -0.044091056, -0.036608716, -0.030612560, -0.025794361, -0.021924841, -0.018225052, -0.015247852, -0.012386003]
    
    TAU_out=np.interp(water_vapour, H2O_list, TAU_list)
    P0_out =np.interp(water_vapour, H2O_list, P0_list)
    P1_out =np.interp(water_vapour, H2O_list, P1_list)
    P2_out =np.interp(water_vapour, H2O_list, P2_list)
    
    vangle=np.radians(vangle)
    atmos_trans=np.exp(-TAU_out/np.cos(np.radians(P0_out + vangle*P1_out + vangle**2*P2_out)))
    return atmos_trans

def get_Atmostrans_MODIS(vangle, water_vapour):
    """
    Look Up Table to calculate the atmospheric transmissivity for MODIS, using
    an user inputted angle, and water vapour.

    Parameters
    ----------
    vangle : float
        Satellite View Angle.
    water_vapour : float
        Water vapour in kg /m2 units.

    Returns
    -------
    atmos_trans : float
        Atmospheric Transmissivity.
    """
    H2O_list=[5, 10, 20, 30, 40, 50, 60]
    TAU_list=[0.1381, 0.1399, 0.1445, 0.1504, 0.1576, 0.166, 0.1757]
    P0_list=[-0.0651622150, -0.0643401720, -0.0629781890, -0.0615412040, -0.0579979620, -0.0603463030, -0.0671102320]
    P1_list=[+1.1533658380, +1.1517193010, +1.1485233610, +1.1450919980, +1.1382583090, +1.1407634180, +1.1498490810]
    P2_list=[-0.0739194221, -0.0732436675, -0.0718538574, -0.0703834389, -0.0675955496, -0.0683776910, -0.0717979974]
    
    TAU_out=np.interp(water_vapour, H2O_list, TAU_list)
    P0_out =np.interp(water_vapour, H2O_list, P0_list)
    P1_out =np.interp(water_vapour, H2O_list, P1_list)
    P2_out =np.interp(water_vapour, H2O_list, P2_list)
    
    vangle=np.radians(vangle)
    atmos_trans=np.exp(-TAU_out/np.cos(np.radians(P0_out + vangle*P1_out + vangle**2*P2_out)))
    return atmos_trans

########## Date Functions
def doy_finder(yr, mn, dy):
    """Calculates the day of year from user input year, month,day.
    Input: yr: year (yyyy)
    Input: mn: month (mm)
    Input: dy: day (dd)

    Output: doy: Day of Year (1st Jan = 1)"""

    doy = datetime.date(int(yr), int(mn), int(dy)).timetuple().tm_yday
    return doy

def date_finder(yr, doy):
    """
    Calculates the month and day from user input of year and day of year.

    Parameters
    ----------
    yr : String
        Year of Acquisition.
    doy : Integar
        Day of Year of Acquistion.

    Returns
    -------
    M_NUM : String
        Two-digit Month.
    DAY : String
        Two-digit Day.

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

########### MODIS Tile Location ###########
def lat_lon_to_modis(lat, lon):
    """
    Function that calculates the horizontal and vertical tile numbers of a 
    given latitude and longitude pairing.

    Parameters
    ----------
    lat : float
        Latitude.
    lon : float
        Longitude.

    Returns
    -------
    h, v : integar
        Horizontal and Vertical tile number.

    """
    p_modis_grid = Proj('+proj=sinu +R=6371007.181 +nadgrids=@null +wktext')
    VERTICAL_TILES = 18
    HORIZONTAL_TILES = 36
    EARTH_RADIUS = 6371007.181 # Radius of the Earth in KM
    EARTH_WIDTH = 2 * np.pi * EARTH_RADIUS
    
    TILE_WIDTH = EARTH_WIDTH / HORIZONTAL_TILES
    TILE_HEIGHT = TILE_WIDTH
    
    x, y = p_modis_grid(lon, lat)
    h = (EARTH_WIDTH * .5 + x) / TILE_WIDTH
    v = -(EARTH_WIDTH * .25 + y - (VERTICAL_TILES - 0) * TILE_HEIGHT) / TILE_HEIGHT
    return int(h), int(v)

def tile_def(lat, lon):
    """
    Function to get the MODIS Tile for a given Latitude and Longitude.

    Parameters
    ----------
    lat : float
        Latitude.
    lon : float
        Longitude.

    Returns
    -------
    tile_str : string
        String of MODIS Sin Tiles, in "hXXvYY" format, where XX is the 
        horizontal tile number, and YY is the vertical tile number.
    """
    h, v = lat_lon_to_modis(lat, lon)
    if h < 10:
        h_str = 'h0'+str(h)
    else:
        h_str = 'h'+str(h)
    if v < 10:
        v_str = 'v0'+str(v)
    else:
        v_str = 'v'+str(v)
    
    tile_str = h_str+v_str
    return tile_str
