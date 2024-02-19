# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 10:46:13 2024

@author: k2262276
"""
from pyhdf.SD import SD, SDC

#%%
def extract_hdf_2D(filename, var_name):
    """ Extracts the 2d Array data from the hdf4 file """
    hdf = SD(filename, SDC.READ)
    var = hdf.select(var_name)[:,:]
    return var

def create_latlon_hdf(old_hdf, new_hdf):
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
    lat = extract_hdf_2D(old_hdf, 'Latitude')
    lon = extract_hdf_2D(old_hdf, 'Longitude')
    
    # Saves attribute "title" to be the filename of the MXD07 file used
    att = filehdf.attr('title')
    att.set(SDC.CHAR, f"Taken from {new_hdf.split('99.')}")
    
    # Create the latitude dataset : set its name, data type and dimensions
    sds_name = "Latitude"
    sds = filehdf.create(sds_name, SDC.FLOAT32, lat.shape)
    sds.dim(0).setname("Cell_Along_Swath:mod07")
    sds.dim(1).setname("Cell_Across_Swath:mod07")
    sds[:] = lat    # Write the latitude to it
    sds.endaccess() # Terminate access to the data set
    
    # Create the longitude dataset : set its name, data type and dimensions
    sds_name = "Longitude"
    sds = filehdf.create(sds_name, SDC.FLOAT32, lon.shape)
    sds.dim(0).setname("Cell_Along_Swath:mod07")
    sds.dim(1).setname("Cell_Across_Swath:mod07")
    sds[:] = lon    # Write the longitude to it
    sds.endaccess() # Terminate access to the data set
    
    # Close the file
    filehdf.end()

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

