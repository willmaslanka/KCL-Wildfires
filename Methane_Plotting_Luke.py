# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 16:16:22 2025

Example of loading in the S5P CH4 product, and plotting

@author: WMaslanka
"""
from area import area
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cartopy.feature as cf
import cartopy.crs as ccrs

"""
You won't have this library, as it's one that I've made up of 
useful functions that I always use. It's not vital for this script, however.
"""
import CDF 

#%% Script Specific Functions
def S5P_pixel_size(S5P_lat_bnds, S5P_lon_bnds):
    """
    Calculated the pixel size of the Sentinel-5p data, assuming the coordinates
    are in WGS84 format (Sentinel-5P documenation). Loads in the S5P pixel 
    corners, generates a polygon, and calculates the area in m2.

    Parameters
    ----------
    S5P_lat_bnds : Array
        Array of pixel corners (latitudes)
    S5P_lon_bnds : Array
        Array of pixel corners (longitude)

    Returns
    -------
    pixel_m2 : Array
        Array of pixel areas in m2.

    """
    # Create Empty Array
    pixel_m2 = np.zeros((S5P_lon_bnds.shape[0], S5P_lon_bnds.shape[1]))
    
    # Looping through the pixels
    for y_idx in range(S5P_lon_bnds.shape[1]): 
        for x_idx in range(S5P_lon_bnds.shape[0]):
            
            # Extract Corner Lat/Lons
            pnt_lat = S5P_lat_bnds[x_idx,y_idx,:]
            pnt_lon = S5P_lon_bnds[x_idx,y_idx,:]
            pixel_points = []
            for i in range(pnt_lat.shape[0]):
                crnr_points = []
                crnr_points.append(pnt_lon[i])
                crnr_points.append(pnt_lat[i])
                pixel_points.append(crnr_points)
            pixel_points.append(pixel_points[0])
            
            # Make polygon of pixel, and calculate area
            obj = {'type':'Polygon','coordinates':[pixel_points]}
            pixel_m2[x_idx, y_idx] = area(obj)
    
    return pixel_m2

#%% Directories and Filenames
""" This is just how I found/loaded in the files - you should be able to just 
put the variable 'fname' to be equal to the absolute filepath of the SP5 file
"""

data_dir = 'C:/Users/k2262276/Documents/Saddleworth2018/Sentinel-5P/'
fnames = CDF.glob_file_checker(data_dir, '*.nc')
fname = fnames[1]


layer_num = -1 # This is for defining which layer in the S5P file I want.
# I believe it is in descending order (layer 0 is the topmost layer, layer -1 is the bottom)

#%% Extracting data from the S5P file
"""
The snippet of code below is what I've used to extract the data from the S5P file.
Basically, the variables are either in either:
    (time, scanline, ground_pixel, layer) shape [A]
    (time, scanline, ground_pixel) shape [B]
    (time, scanline, ground_pixel, corner) shape [C]
    
Regardless, 'time' is always 1 (as it is only scanned once per swath), so the 
np.squeeze function removes that dimension.
"""

ds = nc.Dataset(fname)
ds_prod = ds['PRODUCT']
lat = np.squeeze(ds_prod['latitude'][:]) # [B]
lon = np.squeeze(ds_prod['longitude'][:]) # [B]
mth = np.squeeze(ds_prod['methane_mixing_ratio_bias_corrected'][:]) # [B]
lat_bnds = np.squeeze(ds_prod['SUPPORT_DATA']['GEOLOCATIONS']['latitude_bounds'][:]) # [C]
lon_bnds = np.squeeze(ds_prod['SUPPORT_DATA']['GEOLOCATIONS']['longitude_bounds'][:]) # [C]
alt = np.squeeze(ds_prod['SUPPORT_DATA']['INPUT_DATA']['altitude_levels'][:]) # [A]
mpa = np.squeeze(ds_prod['SUPPORT_DATA']['INPUT_DATA']['methane_profile_apriori'][:]) # [A]

#%% Conversion Factors
# This is my script to calculate the pixel size from the lat/lon corners.
pixel_m2 = S5P_pixel_size(lat_bnds, lon_bnds)

# This is the conversion factor to convert from mol m-2 to g
conv_fact = pixel_m2 * (12 + 1*4) # CH4 = (12 + (1*4)) g/mol

#%% Plotting
#%%% Plotting the layer altitude
fig, axs = plt.subplots(1,2, figsize=(15,8), subplot_kw={'projection': ccrs.PlateCarree()})
axs[0].set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
pcm = axs[0].pcolormesh(lon, lat, alt[:,:,layer_num], transform=ccrs.PlateCarree())
axs[0].add_feature(cf.BORDERS, linestyle=':')
axs[0].add_feature(cf.COASTLINE, linestyle='-')

axs[1].set_extent([-10, 5, 48, 63], ccrs.PlateCarree())
pcm = axs[1].pcolormesh(lon, lat, alt[:,:,layer_num], transform=ccrs.PlateCarree())
axs[1].add_feature(cf.BORDERS, linestyle=':')
axs[1].add_feature(cf.COASTLINE, linestyle='-')

plt.colorbar(pcm, ax=axs, orientation='horizontal', 
             label=f"Layer {layer_num} Height above Geiod [m]")
plt.suptitle(fname.split('/')[-1])

#%%% Plotting the Methane Mixing Ratio
fig, axs = plt.subplots(1,2, figsize=(15,8), subplot_kw={'projection': ccrs.PlateCarree()})
axs[0].set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
pcm = axs[0].pcolormesh(lon, lat, mth, transform=ccrs.PlateCarree())
axs[0].add_feature(cf.BORDERS, linestyle=':')
axs[0].add_feature(cf.COASTLINE, linestyle='-')

axs[1].set_extent([-10, 5, 48, 63], ccrs.PlateCarree())
pcm = axs[1].pcolormesh(lon, lat, mth, transform=ccrs.PlateCarree())
axs[1].add_feature(cf.BORDERS, linestyle=':')
axs[1].add_feature(cf.COASTLINE, linestyle='-')

plt.colorbar(pcm, ax=axs, orientation='horizontal',
             label=r"CH4 mixing ratio (bias corrected) [x10$^{-9}$]")
plt.suptitle(fname.split('/')[-1])

#%%% Plotting the Methane Profile Apriori
fig, axs = plt.subplots(1,2, figsize=(15,8), subplot_kw={'projection': ccrs.PlateCarree()})
axs[0].set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
pcm = axs[0].pcolormesh(lon, lat, mpa[:,:,layer_num], transform=ccrs.PlateCarree())
axs[0].add_feature(cf.BORDERS, linestyle=':')
axs[0].add_feature(cf.COASTLINE, linestyle='-')

axs[1].set_extent([-10, 5, 48, 63], ccrs.PlateCarree())
pcm = axs[1].pcolormesh(lon, lat, mpa[:,:,layer_num], transform=ccrs.PlateCarree())
axs[1].add_feature(cf.BORDERS, linestyle=':')
axs[1].add_feature(cf.COASTLINE, linestyle='-')

plt.colorbar(pcm, ax=axs, orientation='horizontal',
             label=fr"Mole Content of Methane in Layer {layer_num} [mol m$^{-2}$]")
plt.suptitle(fname.split('/')[-1])

#%%% Plotting the Methane Mass for given layer
fig, axs = plt.subplots(1,2, figsize=(15,8), subplot_kw={'projection': ccrs.PlateCarree()})
axs[0].set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
pcm = axs[0].pcolormesh(lon, lat, mpa[:,:,layer_num]*conv_fact, transform=ccrs.PlateCarree())
axs[0].add_feature(cf.BORDERS, linestyle=':')
axs[0].add_feature(cf.COASTLINE, linestyle='-')

axs[1].set_extent([-10, 5, 48, 63], ccrs.PlateCarree())
pcm = axs[1].pcolormesh(lon, lat, mpa[:,:,layer_num]*conv_fact, transform=ccrs.PlateCarree())
axs[1].add_feature(cf.BORDERS, linestyle=':')
axs[1].add_feature(cf.COASTLINE, linestyle='-')

plt.colorbar(pcm, ax=axs, orientation='horizontal',
             label=fr"Mass Methane in Layer {layer_num} [g]")
plt.suptitle(fname.split('/')[-1])



