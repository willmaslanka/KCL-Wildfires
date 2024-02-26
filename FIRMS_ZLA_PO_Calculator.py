# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 09:46:59 2024

@author: k2262276
"""
import pandas as pd
import numpy as np

def LookAngle_FIRMS(filt_df, sat):
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