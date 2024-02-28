# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 10:46:13 2024

@author: k2262276
"""
from pyhdf.SD import SD, SDC
import numpy as np

#%%
def extract_hdf_2D(filename, var_name):
    """ Extracts the 2d Array data from the hdf4 file """
    hdf = SD(filename, SDC.READ)
    var = hdf.select(var_name)[:,:]
    return var

def create_latlon_hdf(old_latlon_file, fire_mask_file, new_hdf):
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
    fm = extract_hdf_2D(fire_mask_file, 'fire mask') - 3 # Removes Not Processed Flags
    fm[fm == 0] = 64 # Mimics MOD14CRS pixel count
    fm[fm == 1] = 32 # Mimics MOD14CRS pixel count
    fm[fm == 2] = 0 # Mimics MOD14CRS pixel count
    cfm = coarsen_array(fm, 5) # Takes the average pixel count of each 5x5 window
    
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
    sds[:] = lon    # Write the longitude to it
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
old_fname = 'filepath-to-MOD06'
product = 'MODO6_L2'
proxy_save_dir = 'directory-to-save-into'
doy = 'day-of-year-value'

#%% Proxy Product Naming
if product[1].upper() == 'O':
    proxy_product = 'MOD99'
else:
    proxy_product = 'MYD99'
    
#%% Script
f_id = old_fname.split(f'{product}.')[-1]
new_hdf_fname = f"{proxy_save_dir}{doy}/{proxy_product}.{f_id}"

# Create daily MXD99 file
create_latlon_hdf(old_fname, new_hdf_fname)

