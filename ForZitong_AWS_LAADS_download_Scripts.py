# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 12:09:35 2024

@author: k2262276
"""

#%% Libraries
import requests, subprocess, os, datetime, sys
from bs4 import BeautifulSoup # Library to scrap information from online
from tqdm import tqdm # Library to create progress bars for loops/functions
import s3fs # Module to interface with Amazon Simple Storage Service (S3)
from packaging.version import parse # Library to access core utilities for Python packages

#%% UDF
def ladsweb_download(collection, product, year, month, day, 
                     start_time, end_time, save_dir, token, ext='nc'):
    """
    Function to download data directly from the LAADS DAAC from NASA.
    Requires an earth data token to download.
    Steps for token generation can be found here:
        https://urs.earthdata.nasa.gov/documentation/for_users/user_token
        
    Will download the files to a user defined directory.

    Parameters
    ----------
    collection : String
        Data Collection Bucket. Collection varies depending on what product / 
        instrument is required. 
        L1B Swath Calibrated Radiences can be found in Collection:
            '5200' for SNPP (VNP) and NOAA-21/JPSS2 (VJ2)
            '5201' for NOAA-20/JPSS1 (VJ1)
        URL of all collections can be found here:
            https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/
    product : String
        Product Abbrieviation.
        The following are the appropriate Multispectral and geolocation products:
                VNP02MOD - VIIRS/NPP Moderate Resolution 6-Min L1B Swath 750m
                VNP02IMG - VIIRS/NPP Imagery Resolution 6-Min L1B Swath 375m
                VNP03MOD - VIIRS/NPP Moderate Resolution Terrain-Corrected Geolocation 6-Min L1 Swath 750m
                VNP03IMG - VIIRS/NPP Imagery Resolution Terrain-Corrected Geolocation 6-Min L1 Swath 375m
                [Data Starts on Day Of Year 019, 2012]
                
                VJ102MOD - VIIRS/JPSS1 Moderate Resolution 6-Min L1B Swath 750m
                VJ102IMG - VIIRS/JPSS1 Imagery Resolution 6-Min L1B Swath 375m
                VJ103MOD - VIIRS/JPSS1 Moderate Resolution Terrain-Corrected Geolocation 6-Min L1 Swath 750m
                VJ103IMG - VIIRS/JPSS1 Imagery Resolution Terrain-Corrected Geolocation 6-Min L1 Swath 375m
                [Data Starts on Day of Year 005, 2018]
                
                VJ202MOD - VIIRS/JPSS2 Moderate Resolution 6-Min L1B Swath 750m
                VJ202IMG - VIIRS/JPSS2 Imagery Resolution 6-Min L1B Swath 375m
                VJ203MOD - VIIRS/JPSS2 Moderate Resolution Terrain-Corrected Geolocation 6-Min L1 Swath 750m
                VJ203IMG - VIIRS/JPSS2 Imagery Resolution Terrain-Corrected Geolocation 6-Min L1 Swath 375m
                [Data Starts on Day of Year 041, 2023]

    year : String
        String of year required (YYYY).
    month: String
        String of month required (MM)
    day: String
        String of day required (DD)
    start_time : String
        Start period, in four digits (HHMM).
    end_time : String
        End period, in four digits (HHMM).
    save_dir : String
        Directory to save the data to.
    token : String
        Earthdata Token.
    ext : String
        Extension needed. Default is 'nc'.
    """
    
    def doy_finder(yr, mn, dy):
        """Calculates the day of year from user input year, month,day.
        Input: yr: year (yyyy)
        Input: mn: month (mm)
        Input: dy: day (dd)

        Output: doy: Day of Year (1st Jan = 1)"""

        doy = datetime.date(int(yr), int(mn), int(dy)).timetuple().tm_yday
        
        if doy > 99:
            s_doy = str(doy)
        elif doy > 9:
            s_doy = '0'+str(doy)
        else:
            s_doy = '00'+str(doy)
        return s_doy

    
    def listFD(url, ext=''):
        # scraps page for files
        page = requests.get(url).text
        soup = BeautifulSoup(page, 'html.parser')
        return [url + '/' + node.get('href') 
                for node in soup.find_all('a') 
                if node.get('href').endswith(ext)]
    
    def runcmd(cmd, verbose = False, *args, **kwargs):
        # Script for running command line functions from python
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
    
    doy = doy_finder(year, month, day)
    # Creating search URL
    master_url = 'http://ladsweb.modaps.eosdis.nasa.gov/archive/allData/'
    search_url = f'{master_url}{collection}/{product}/{year}/{doy}/'

    # Find all files within search times
    files = []
    for file in listFD(search_url, ext):
        #for swath in file:
        acq_time = file.split('///')[-1].split('.')[2]
        if (int(acq_time) >= int(start_time)) and (int(acq_time) <= int(end_time)):
               files.append(file)
    
    files = list(set(files)) # Removing duplicate files
    
    # looping
    for file in tqdm(files, unit='files', 
                     bar_format="{desc}Downloading:{percentage:3.0f}%|{bar}|{n_fmt}/{total_fmt} [{elapsed}<{remaining}]"):
        url = file.split('///')[0] + '/' + file.split('///')[-1]
        if not os.path.isfile(save_dir + url.split('/')[-1]):
            # If the file does not exist in the directory, download it
            runcmd(f'wget -e robots=off -m -np -q -nd -R .html,.tmp -nH --cut-dirs=3 {url} --header "Authorization: Bearer {token}" -P {save_dir}', verbose = True)    

def AWS_VIIRS_download(sat, product, year, month, day, 
                       start_time, end_time, save_dir):
    """
    

    Parameters
    ----------
    sat : String
        Satellite required. Currently, this only works for S-NPP and NOAA-20.
            For S-NPP, the string is 'SNPP'
            For NOAA-20/JPSS1, the string is 'NOAA20'
    product : String
        The Product that you are after.
            Options are either 'AOD' or 'ADP'
            'AOD': Aerosol Optical Depth
            'ADP': Aerosol Detection
    year : String
        String of year required (YYYY).
    month: String
        String of month required (MM)
    day: String
        String of day required (DD)
    start_time : String
        Start period, in four digits (HHMM).
    end_time : String
        End period, in four digits (HHMM).
    save_dir : String
        Directory to save the data to.
    """
    
    def get_product_path(product, processing, sat):
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
        product_path = sat + '_VIIRS_' + product_name + suffix
        return product_path
    
    def aws_viirs_list(year, month, day, start_hour, start_min, end_hour, end_min, satellite, product, processing):
        """ Function to get list of VIIRS data from AWS server for the given time period.
        yyyy/mm/dd start time (hh:mm) to end time (hh:mm)
        for a given satellite, product, and processing level."""
        
        # Construct aerosol product path for AWS NODD
        product_path = get_product_path(product, processing, satellite)

        # Access AWS using anonymous credentials
        aws = s3fs.S3FileSystem(anon=True)
        
        # Create list of file names for entire day (~550)
        day_files = aws.ls('noaa-jpss/' + satellite + '/VIIRS/' + product_path + '/' + year + '/' + month + '/' + day + '/', refresh=True)  
        
        # Create list of subsetted file names that fall within specified time period(s)
        data = []
        start_datetime = datetime.datetime(int(year), int(month), int(day), int(start_hour), int(start_min))
        end_datetime = datetime.datetime(int(year), int(month), int(day), int(end_hour), int(end_min))
        for file in day_files:
            file_daytime = datetime.datetime.strptime(file.split('_')[-3][1:13], '%Y%m%d%H%M')
            if file_daytime >= start_datetime and file_daytime <= end_datetime:
                data.append(file)
        return data
    
    if datetime.datetime(int(year), int(month), int(day)) <= datetime.datetime(2022, 12, 31):
        processing = 'reprocessed'
    else:
        processing = 'operational'
    
    data = aws_viirs_list(year, month, day, 
                          start_time[:2], start_time[2:], 
                          end_time[:2], end_time[2:], 
                          sat, product, processing)
    print('data found')
    
    if len(data) > 0:   
        # Access AWS using anonymous credentials
        aws = s3fs.S3FileSystem(anon=True)
        if parse(sys.version.split(' ')[0]) < parse('3.9'):
                sys.stdout.flush()
        for name in tqdm(data, unit='files', bar_format="{desc}Downloading:{percentage:3.0f}%|{bar}|{n_fmt}/{total_fmt} [{elapsed}<{remaining}]"):
            print(name)
            # Set save_path + file_name as pathlib.Path object and convert to string
            #full_path = str(save_path / name.split('/')[-1])
            full_path = save_dir + str(name.split('/')[-1])
            # Download file from AWS archive
            aws.get(name, full_path)
            file_size = aws.size(name)
            print(name.split('/')[-1], ' (', int(file_size/1.0E6), ' MB)', sep='')
        print('Download complete!')
    else:
        print('No files retrieved. Check settings and try again.')

#%% Search Criteria - ladswed_download
collection = '5200'
product = 'VNP14'
year = '2022'
month = '01'
day = '01'
start_time = '0100'
end_time = '0110'
ext = 'nc'
earthdata_token = 'personal-earthdata-token'
save_dir = 'directory-to-save'

# Download Function
ladsweb_download(collection, product, year, month, day,
                 start_time, end_time, save_dir, earthdata_token)

#%% Search Criteria - AWS_VIIRS_download
sat = 'SNPP'
product = 'AOD'
year = '2023'
month = '01'
day = '01'
start_time = '0100'
end_time = '0115'
save_dir = 'directory-to-save'

AWS_VIIRS_download(sat, product, year, month, day,
                   start_time, end_time, save_dir)