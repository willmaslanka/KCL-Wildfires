# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 09:46:59 2024

@author: k2262276
"""
import pandas as pd
import numpy as np
import datetime
import netCDF4 as nc

def LookAngle_FIRMS(filt_df: pd.DataFrame, sat: str) -> pd.DataFrame:
    """
    Parameters
    ----------
    filt_df : DataFrame
        Pandas Dataframe: taken from FIRMS, relavent for polar orbiters. 
    sat : String
        Polar Orbiter Identifier - Only value for MODIS and VIIRS so far.

    Returns
    -------
    new_df : DataFrame
        Pandas Dataframe: Copy of input DataFrame, but with view angle 
        calculated

    """
    
    Re = 6378.137 # earth radius, km
    if sat == 'MODIS':
        h = 705.0 # satellite altitude, km
        px = 1 # pixel resolution, km
    elif sat == 'VIIRS':
        px = 0.375 # pixel resolution, km
        h = 833.0 # satellite altitude, km
    r = Re + h # satellite altitude from centre of earth, km
    s = px/h # pixel dimension / altitude
    SLA = [] # Empty Array of Satallite Look Angle

    if sat.upper() == 'MODIS': # If MODIS
        alpha = filt_df.track.values / (r * s)
        beta = filt_df.scan.values / (Re * s)
        for coeff in zip(alpha, beta):
            if coeff[0] * (coeff[1] + 1) / (coeff[1]) <= 1:
                sla = np.degrees(np.arccos(coeff[0] * (coeff[1] + 1) / (coeff[1])))
            else:
                sla = np.degrees(np.arccos(1))
            SLA.append(sla)
        
    else: # If VIIRS
        A = (r * s)
        B = ((filt_df.track.values/A) ** 2)
        C = ((Re/r) ** 2)
        D = (2*filt_df.track.values)
        for coeff in zip(B, D):
            if A * (1 + coeff[0] - C)/coeff[1] <= 1:
                sla = np.degrees(np.arccos(A * (1 + coeff[0] - C)/coeff[1]))
            else:
                sla = np.degrees(np.arccos(1))
            SLA.append(sla)
 
    new_df = filt_df.copy()
    new_df['zla'] = SLA
    
    return new_df

def closest_value(input_list: list, input_value: float) -> float | int:
    """ Gets nearest value in list to input value """
    arr = np.asarray(input_list)
    i = (np.abs(arr - input_value)).argmin()
    return arr[i], i

def point_TCWV_ECMWF(lat_pnt: float,
                     lon_pnt: float,
                     date_pnt: 'datetime',
                     ecmwf_lat: float,
                     ecmwf_lon: float,
                     ecmwf_tme: float,
                     tcwv_pro: float) -> float:
    """ Extracts the TCWV for a given latitude/longitude 
    lat_pnt = latitude at point
    lon_pnt = longitude at point
    date_pnt = datetime at point
    ecmwf_lat = latitude from ECMWF TCWV file,
    ecmwf_lon = longitude from ECMWF TCWV file,
    ecmwf_tme = time from ECMWF TCWV file,
    tcwv_pro = processes TCWV from ECMWF TCWV file
    """
    __, lon_idx = closest_value(ecmwf_lon, lon_pnt)
    __, lat_idx = closest_value(ecmwf_lat, lat_pnt)
    
    atmos_epoch = datetime.datetime(1900,1,1,0,0)
    time_diff = date_pnt - atmos_epoch
    hours_past = time_diff.total_seconds()/60/60
    pnts = int(np.round(hours_past % 1,1)*10)
    
    if pnts == 0:
        tme_idt = np.where(ecmwf_tme == np.floor(hours_past))[0][0]
        
        tcwv_pnt = tcwv_pro[tme_idt, lat_idx, lon_idx]
    else:
        tme_idt1 = np.where(ecmwf_tme == np.floor(hours_past))[0][0]
        tme_idt2 = np.where(ecmwf_tme == np.ceil(hours_past))[0][0]
        
        tcwv_idt1 = tcwv_pro[tme_idt1, lat_idx, lon_idx]
        tcwv_idt2 = tcwv_pro[tme_idt2, lat_idx, lon_idx]

        tcwv_linspace = np.linspace(tcwv_idt1, tcwv_idt2, 11)
        tcwv_pnt = tcwv_linspace[pnts]
    
    return tcwv_pnt

def get_Atmostrans_VIIRS(vangle: float, water_vapour: float) -> float:
    """ Applies MODTRAN4 To VIIRS data 
    vange = view angle, water_vapour = ECMWF Total Column Water Vapour
    atmos_trans = atmospheric transmissivity """
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

def get_Atmostrans_MODIS(vangle: float, water_vapour: float) -> float:
    """ Applies MODTRAN4 To MODIS data 
    vange = view angle, water_vapour = ECMWF Total Column Water Vapour
    atmos_trans = atmospheric transmissivity """
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

def ApplyAtmosCorrection(pre_mask_df: pd.DataFrame,
                         sat: str,
                         ecmwf_lat: float,
                         ecmwf_lon: float,
                         ecmwf_tme: int,
                         tcwv_pro: float) -> pd.DataFrame:
    """ Master Function for applying the Atmospheric Correction to FIRMS data
    pre_mask_df = FIRMS dataframe,
    sat = satellite (MODIS or VIIRS)
    ecmwf_lat = latitude from ECMWF TCWV file,
    ecmwf_lon = longitude from ECMWF TCWV file,
    ecmwf_tme = time from ECMWF TCWV file,
    tcwv_pro = processes TCWV from ECMWF TCWV file
    """
    
    latitude = pre_mask_df.latitude.values
    longitude = pre_mask_df.longitude.values
    year = pre_mask_df.year.values
    month = pre_mask_df.month.values
    day = pre_mask_df.day.values
    hour = pre_mask_df.hour.values
    minute = pre_mask_df.minute.values
    zla = pre_mask_df.zla.values
    
    TAU = []
    for vars in zip(latitude, longitude, year, month, day, hour, minute, zla):
        frp_dt = datetime.datetime(int(vars[2]), int(vars[3]), int(vars[4]), int(vars[5]), int(vars[6]))
        tcwv = point_TCWV_ECMWF(vars[0], vars[1], frp_dt, ecmwf_lat, ecmwf_lon, tcwv_pro, tcwv_pro)
        
        if sat.upper() == 'VIIRS':
            tau = get_Atmostrans_VIIRS(vars[7], tcwv)
        elif sat.upper() == 'MODIS':
            tau = get_Atmostrans_MODIS(vars[7], tcwv)
        
        TAU.append(tau)
    
    new_df = pre_mask_df.copy()
    new_df['tau'] = TAU
    new_df['atm_frp'] = new_df.frp / new_df.tau
    
    return new_df

def extract_tcwv(tcwv_fname: str) -> float|float|float|float:
    """ Extracting TCWV data from ECMWF netCDF file
    tcwv_fname = string to file"""
    
    ds = nc.Dataset(tcwv_fname)
    lat = ds['latitude'][:]
    lon = ds['longitude'][:]
    tme = ds['time'][:] # hours since 1900-01-01 00:00:00.0
    tcwv_raw = ds['tcwv'][:] # kg m**-2
    tcwv_scl = ds['tcwv'].scale_factor
    tcwv_add = ds['tcwv'].add_offset
    ds.close()
    
    tcwv_pro = (tcwv_raw*tcwv_scl) + tcwv_add
    
    return lat, lon, tme, tcwv_pro

