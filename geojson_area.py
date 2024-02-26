# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:26:43 2024

@author: WMaslanka

Need to install the "area" library (https://pypi.org/project/area/)
"""

from area import area
import geojson

fname = 'filepath-of-geojson'

def plume_area_calc(fname: str) -> float:
    """ Calculates the area of the geojson polygon, in Sq. Meters.
    If wanted in Sq. km, then Sq. Meters / (10**6)"""
    with open(fname) as f:
        gj = geojson.load(f)
    return area(gj['features'][0]['geometry'])

area_sq_m = plume_area_calc(fname)
