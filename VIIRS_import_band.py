# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 11:37:57 2023
@author: Will Maslanka

List of Functions used to import VIIRS data into Python. Makes use of the 
netCDF4 and Numpy Packages.

Main Functions used within this python script are the following:
    *** VIIRS_M_Band and VIIRS_I_Band ***

In Theory, If you have the absolute filepath of the VIIRS Moderiate-Resolution 
data, along with the absolute filepaths associated VIIRS Image-Resolution data,
 and VIIRS Image-Resolution Geolocation data, the follwoing functions will give
 you arrays of chosen bands at image-resolution, as well as the geolocation 
information at image-resolution.
"""

#%% Modules
import netCDF4 as nc
import numpy as np

#%% User-Defined-Functions
def VIIRS_M_Band(mref_file, scle='scaled', Red='M05', Green='M04', Blue='M03'):
    """
    Parameters
    ----------
    mref_file : String
        Moderate-resolution Calibrated Radiances.
    scle : String
        Toggle for scaling array ('scaled' allows for scaling via percentage, 
                                  otherwise it is MinMax).
    Red : String
        Band name for what you want in the Red slot of a Colour Composite
        Default value: 'M05'
    Green : String
        Band name for what you want in the Green slot of a Colour Composite
        Default value: 'M04'
    Blue : String
        Band name for what you want in the Blue slot of a Colour Composite
        Default value: 'M03'

    Returns
    -------
    B : Array
        Array of the Reflectances in the Blue Colour Composite Position,
        at Image-Resolution.
    G : Array
        Array of the Reflectances in the Green Colour Composite Position,
        at Image-Resolution.
    R : Array
        Array of the Reflectances in the Red Colour Composite Position,
        at Image-Resolution.

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
    contains = np.in1d(scaled_band,flags).reshape((scaled_band.shape[0],
                                                   scaled_band.shape[1]))
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
    x : Array
        Array to be scaled.
    pcnt : Integar
        How much to trim from the array.

    Returns
    -------
    scle : Array
        Scaled Array.

    """
    scle = (x - np.nanpercentile(x, pcnt)) / (np.nanpercentile(x, 100 - pcnt) 
                                              - np.nanpercentile(x, pcnt))
    
    scle[scle >= 1] = 1
    scle[scle <= 0] = 0
    
    return scle