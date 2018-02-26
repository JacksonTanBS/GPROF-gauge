#!/usr/bin/python

import numpy as np

def haversine(lat1, lon1, lat2, lon2, R = 6372.8):
    """The Haversine function to calculate distance between two sets of 
       lat/lon coordinates. Inputs are in degrees; output is in km."""

    dLat = np.deg2rad(lat2 - lat1)
    dLon = np.deg2rad(lon2 - lon1)
    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)
 
    a = np.sin(dLat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dLon / 2) ** 2
    c = 2 * np.ma.arcsin(np.ma.sqrt(a))
 
    return R * c
