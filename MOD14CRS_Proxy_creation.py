# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 10:46:13 2024

@author: k2262276
"""
from pyhdf.SD import SD, SDC
import numpy as np
import os
from glob import glob

#%%
def extract_hdf_2D(filename, var_name):
    """ Extracts the 2d Array data from the hdf4 file """
    hdf = SD(filename, SDC.READ)
    var = hdf.select(var_name)[:,:]
    return var

def glob_file_checker(directory, search_string):
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

def create_MOD14CRS_hdf(old_latlon_file, fire_mask_file, new_hdf):
    """ 
    Writes the new hdf4 file
    
    Specifically, it takes the MXD07 file, and reads the latitude and longitude 
    variables from it. It then creates a new file (proxy product name of MXD99) 
    that contains just the latitude and longitude data. The file name after the 
    product name is taken from the parent MXD07 file. The MXD99 file also has a 
    title, stating which MXD07 file it was taken from.
    """
    # Creates new empty  hdf4 file
    filehdf  = SD(new_hdf, SDC.WRITE | SDC.CREATE)
    
    # Extracts the latitude and longitude from the MXD07 file
    lat = extract_hdf_2D(old_latlon_file, 'Latitude')
    lon = extract_hdf_2D(old_latlon_file, 'Longitude')
    
    # Extracts and sorts out the 
    if fire_mask_file.upper != 'NOTEXIST':
        fm = extract_hdf_2D(fire_mask_file, 'fire mask') - 3 # Removes Not Processed Flags
        fm[fm == 0] = 64 # Mimics MOD14CRS pixel count
        fm[fm == 1] = 32 # Mimics MOD14CRS pixel count
        fm[fm == 2] = 0 # Mimics MOD14CRS pixel count
        cfm = coarsen_array(fm, 5) # Takes the average pixel count of each 5x5 window
    else:
        cfm = np.zeros((lat.shape[0], lat.shape[1]))
        
    cfm = cfm[:,:-1].astype('uint8')
    
    # Saves attribute "title" to be the filename of the MXD07 file used
    att = filehdf.attr('title')
    att.set(SDC.CHAR, f"Taken from {new_hdf.split('99.')}")
    
    # Create the latitude dataset : set its name, data type and dimensions
    sds_name = "Latitude"
    sds = filehdf.create(sds_name, SDC.FLOAT32, lat.shape)
    sds.dim(0).setname("number_of_coarse_resolution_scan_lines")
    sds.dim(1).setname("coarse_resolution_pixels_per_scan_line")
    sds[:] = lat    # Write the latitude to it
    sds.endaccess() # Terminate access to the data set
    
    # Create the longitude dataset : set its name, data type and dimensions
    sds_name = "Longitude"
    sds = filehdf.create(sds_name, SDC.FLOAT32, lon.shape)
    sds.dim(0).setname("number_of_coarse_resolution_scan_lines")
    sds.dim(1).setname("coarse_resolution_pixels_per_scan_line")
    sds[:] = lon    # Write the longitude to it
    sds.endaccess() # Terminate access to the data set
    
    # Create the longitude dataset : set its name, data type and dimensions
    sds_name = "CoarseFireCounts"
    sds = filehdf.create(sds_name, SDC.UINT8, cfm.shape)
    sds.dim(0).setname("number_of_coarse_resolution_scan_lines")
    sds.dim(1).setname("coarse_resolution_pixels_per_scan_line")
    sds[:] = cfm    # Write the longitude to it
    sds.endaccess() # Terminate access to the data set
    
    # Close the file
    filehdf.end()
    
def coarsen_array(array, coarseness=5):
    """ 
    Coarsens the array by a factor of coarseness.
    """
    shape = np.array(array.shape, dtype=int)
    new_shape = coarseness * np.ceil(shape / coarseness).astype(int)

    # Create the zero-padded array and assign it with the old density
    zp_pop_density = np.zeros(new_shape)
    zp_pop_density[:shape[0], :shape[1]] = array

    # Now use the same method as before
    temp = zp_pop_density.reshape((new_shape[0] // coarseness, coarseness,
                                   new_shape[1] // coarseness, coarseness))
    coarse_array = np.round(np.mean(temp, axis=(1,3)) / 32) * 32 # Rounds to nearest 32
    
    return coarse_array

#%% User Variables
#latlon_5m_file = 'filepath-to-MOD06'
#firemask_file = 'filepath-to-MOD14'

#product = 'MODO6_L2'
#proxy_save_dir = 'directory-to-save-into'
#doy = 'day-of-year-value'

latlon_5m_dir = 'C:/Users/k2262276/Documents/DataDownload/CRS_Comparison/MOD07_L2/'
firemask_dir = 'C:/Users/k2262276/Documents/DataDownload/CRS_Comparison/MOD14/'

product = 'MOD07_L2'
proxy_save_dir = 'C:/Users/k2262276/Documents/DataDownload/CRS_Comparison/MOD99/'
doy = '001'

#%% Proxy Product Naming
if product[1].upper() == 'O':
    proxy_product = 'MOD99'
else:
    proxy_product = 'MYD99'
    
#%% Finding all MOD07 files (lat/lon)
MOD07_file = glob_file_checker(f'{latlon_5m_dir}{doy}/', '*')


for latlon_5m_file in MOD07_file: # For every lat/lon file
    time = latlon_5m_file.split(f'{doy}.')[-1].split('.')[0] # Extract Time
    firemask_file = glob_file_checker(f'{firemask_dir}{doy}/',f'*.{time}.*')[0] # Find MOD14 file at same time
    if not firemask_file: # If it doesn't exist, flag
        firemask_file = 'NotExist'
    
    f_id = latlon_5m_file.split(f'{product}.')[-1] # Get Product Extension
    new_hdf_fname = f"{proxy_save_dir}{doy}/{proxy_product}.{f_id}" # Make New Filename for MXD99
    create_MOD14CRS_hdf(latlon_5m_file, firemask_file, new_hdf_fname) # Create MXD99 file