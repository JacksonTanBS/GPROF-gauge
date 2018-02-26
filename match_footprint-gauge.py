#!/usr/bin/python

'''
Matching GPROF/DPR footprints to gauges
Version 0.4.1

This script matches GPROF/DPR footprints to gauge locations. For each month, it 
will select all the GPROF/DPR files for a specified platform and sensor, compare 
them to gauges of a specified netowrk, and record the relevant information in a 
text file. Each match contains the dates, orbit number, the indices of the GPROF 
pixel, the gauges matching the footprint, and the distances between the gauges
and the footprint.

Jackson Tan (jackson.tan@nasa.gov)
2017 06 21
'''


import numpy as np
import sys
import os
import h5py
import tarfile
from calendar import monthrange
from datetime import datetime
from tools import haversine


#--- Preliminaries ---#


year, month = int(sys.argv[1]), int(sys.argv[2])
platform, sensor, network = sys.argv[3], sys.argv[4], sys.argv[5]
datapath = '/home/jackson/Work/Data/footprint-gauge/'
gaugepath = '/media/jackson/Vault/'
if sensor == 'Ku':
    satpath = '/media/jackson/Vault/DPR/'
else:
    satpath = '/media/jackson/GPROF/V05/'


#--- Settings ---#


maxnetworkdist = 200  # max distance (km) between network center and footprint
maxgaugedist = 30     # max distance (km) between gauge and footprint
network_center = {'PCMK': (38.075, -75.57), 'WEGN': (46.90, 15.90),
                  'OLPX': (47, -123.5), 'WGEW': (31.71, -110.05)}


#--- Defining the required functions ---#


def getLocalFiles(path, year, month, platform, sensor):
    
    files = []

    for day in range(1, monthrange(year, month)[1] + 1):
        inpath = '%s%4d/%02d/%02d/' % (path, year, month, day)        
        for file in sorted(os.listdir(inpath)):
            if '%s.%s' % (platform, sensor) in file:
                files.append('%s%s' % (inpath, file))

    return files

def readPCMKlatlon(year, month):

    pcmkpath = '%sGAG/PCMK/combined/' % gaugepath
    gauges = list(range(1, 26))
    buckets = ('A', 'B')

    lat_gauges = []
    lon_gauges = []
    gauges_fixed = []

    for gauge in gauges:

        lat_gauge = []
        lon_gauge = []

        with open('%sIOWA-NASA00%02d_A.gag' % (pcmkpath, gauge), 'r') as f:
            data = [line.split() for line in f.readlines()]

        for ll in data[2:]:
            if int(ll[0]) == year and int(ll[1]) == month:
                lat_gauge.append(float(ll[8]))
                lon_gauge.append(float(ll[9]))

        if len(lat_gauge) == 0:
            continue

        if ((max(lat_gauge) - min(lat_gauge) < 0.01) and 
            (max(lon_gauge) - min(lon_gauge) < 0.01)):
            lat_gauges.append(np.mean(lat_gauge))
            lon_gauges.append(np.mean(lon_gauge))
            gauges_fixed.append(gauge)
            
    return lat_gauges, lon_gauges, gauges_fixed

def readWEGNlatlon(year, month):

    wegpath = '%sWegenerNet/Station/' % gaugepath
    gauges = list(range(1, 152))

    lat_gauges = {}
    lon_gauges = {}
    gauges_recorded = []

    with open('%swegenernet_stations.csv' % wegpath, 'r') as f:
        for rawline in f.readlines()[1:]:
            line = rawline.split(',')
            gauge = int(line[0])
            if gauge not in gauges:
                continue
            startyear = int(line[1][0 : 4])
            startmonth = int(line[1][5 : 7])
            if (datetime(year, month, 1) == datetime(startyear, startmonth, 1)
                and gauge in gauges_recorded):
                gauges_recorded.remove(gauge)  # remove if stn moved this month
            elif datetime(year, month, 1) > datetime(startyear, startmonth, 1):
                lon_gauges[gauge] = float(line[3])
                lat_gauges[gauge] = float(line[2])
                gauges_recorded.append(gauge)

    gauges_recorded = list(set(gauges_recorded))
    lat_gauges = [lat_gauges[gauge] for gauge in gauges_recorded]
    lon_gauges = [lon_gauges[gauge] for gauge in gauges_recorded]

    return lat_gauges, lon_gauges, gauges_recorded

def readOLPXlatlon(year, month):

    olpxpath = '%sGAG/OLPX/combined/' % gaugepath
    gauges = [26, 27, 28, 29, 31, 32, 33, 34, 35, 37, 39, 40, 
              41, 42, 43, 44, 45, 101]
    buckets = ('A', 'B')

    lat_gauges = []
    lon_gauges = []
    gauges_fixed = []

    for gauge in gauges:

        lat_gauge = []
        lon_gauge = []

        with open('%sOLYMPEx-NASA%04d_A.gag' % (olpxpath, gauge), 'r') as f:
            data = [line.split() for line in f.readlines()]

        for ll in data[2:]:
            if int(ll[0]) == year and int(ll[1]) == month:
                lat_gauge.append(float(ll[8]))
                lon_gauge.append(float(ll[9]))

        if len(lat_gauge) == 0:
            continue

        if ((max(lat_gauge) - min(lat_gauge) < 0.01) and 
            (max(lon_gauge) - min(lon_gauge) < 0.01)):
            lat_gauges.append(np.mean(lat_gauge))
            lon_gauges.append(np.mean(lon_gauge))
            gauges_fixed.append(gauge)
            
    return lat_gauges, lon_gauges, gauges_fixed

def readWGEWlatlon(year, month):

    from glob import glob

    wgewpath = '%sWGEW/' % gaugepath

    gauges_available = [int(ii.split('_')[1].split('.')[0]) for ii in 
                        glob('%swgew_*.txt' % wgewpath)]

    gauges, lat_gauges, lon_gauges = [], [], []
    with open('%scoordinates.txt' % wgewpath, 'r') as f:
        for line in f.readlines()[1:]:
            l = line.split()
            if int(l[1]) in gauges_available:
                gauges.append(int(l[1]))
                lon_gauges.append(float(l[2]))
                lat_gauges.append(float(l[3]))
            
    return lat_gauges, lon_gauges, gauges

def calcOP(all_files, lat_network, lon_network):

    files = []
    coords = []
    gauges = []
    distances = []

    for file in all_files:

        with h5py.File(file, 'r') as f:
            if sensor == 'Ku':
                lats_sat = f['NS/Latitude'][:]
                lons_sat = f['NS/Longitude'][:]
            else:
                lats_sat = f['S1/Latitude'][:]
                lons_sat = f['S1/Longitude'][:]

        if len(lats_sat) == 0 or len(lons_sat) == 0:
            continue    # skip problematic files

        if np.min(haversine(lats_sat, lons_sat, 
                            lat_network, lon_network)) <= maxnetworkdist:

            # for each gauge, get the indices of the closest satellite pixel
            coords_gauge = {}
            dists_gauge = {}
            for gg, gauge in enumerate(gauges_fixed):
                dists = haversine(lats_sat, lons_sat, 
                                  lat_gauges[gg], lon_gauges[gg])
                dmin = np.ma.where(dists < maxgaugedist)
                coords_gauge[gauge] = [coord for coord in zip(*dmin)]
                for coord in zip(*dmin):    # record gauge-footprint distance
                    dists_gauge[gauge, coord] = dists[coord]

            coords_all = set([jj for ii in coords_gauge.values() for jj in ii])

            # for each footprint, record the file, coords. and gauge no.
            for coord_min in sorted(coords_all):
                files.append(file)
                coords.append(coord_min)
                gauges.append([gauge for gauge, coord in coords_gauge.items() 
                               if coord_min in coord])
                distances.append([dists_gauge[gauge, coord_min] 
                                  for gauge, coord in coords_gauge.items() 
                                  if coord_min in coord])

    return files, coords, gauges, distances

def readGaugeLatLon(network, *args):
    if   network == 'PCMK':
        return readPCMKlatlon(*args)
    elif network == 'WEGN':
        return readWEGNlatlon(*args)
    elif network == 'OLPX':
        return readOLPXlatlon(*args)
    elif network == 'WGEW':
        return readWGEWlatlon(*args)
    else:
        sys.exit('Error: network not specified.')


#--- Perform the matching ---#


all_files = getLocalFiles(satpath, year, month, platform, sensor)
lat_gauges, lon_gauges, gauges_fixed = readGaugeLatLon(network, year, month)
files, coords, gauges, distances = calcOP(all_files, *network_center[network])


#--- Write the output ---#


inpath = '%s%s/' % (datapath, network)
if not os.path.exists(inpath): os.makedirs(inpath)

# extract the dates and orbit numbers
orbits = [file.split('/')[-1].split('.')[5] for file in files]
dates = [file.split('/')[-4 : -1] for file in files]

# write individual files
file1 = '%s.%s_%4d%02d_footprints.txt' % (platform, sensor, year, month)
with open('%s%s' % (inpath, file1), 'w') as f:
    fmt = '{:4s} {:2s} {:2s} {:6s} {:4d} {:4d}\n'
    for dd, oo, cc in zip(dates, orbits, coords):
        f.write(fmt.format(dd[0], dd[1], dd[2], oo, cc[0], cc[1]))
file2 = '%s.%s_%4d%02d_gauges.txt' % (platform, sensor, year, month)
with open('%s%s' % (inpath, file2), 'w') as f:
    for gg in gauges:
        fmt = '{:3d} ' * len(gg) + '\n'
        f.write(fmt.format(*gg))
file3 = '%s.%s_%4d%02d_distances.txt' % (platform, sensor, year, month)
with open('%s%s' % (inpath, file3), 'w') as f:
    for dd in distances:
        fmt = '{:5.2f} ' * len(dd) + '\n'
        f.write(fmt.format(*dd))

# put them into a compressed tarball
filename = '%s.%s_%4d%02d.tgz' % (platform, sensor, year, month)
with tarfile.open('%s%s' % (inpath, filename), 'w:gz') as f:
    for file in (file1, file2, file3):
        f.add('%s%s' % (inpath, file), arcname = file)
        os.remove('%s%s' % (inpath, file))


'''
Version history
0.4.1    170621    Added WGEW network.
0.4.0    170525    Changed code to work with new WEGN data format.
0.3.1    161123    Fixed minor error in reading WEGN position.
0.3.0    161122    Integrated DPR Ku (match_Ku-gauge.py). Added OLPX.
0.2.4    161121    Modified PCMK coordinates reading function to use combined 
                   codes.
0.2.3    161108    Added handling of masked elements for the rare occasion that
                   the GPROF coordinates are masked.
0.2.2    161107    Added dates to the footprint output to make the reassembly
                   of filename from orbit numbers easier.
0.2.1    161103    Added distance as another output. Also changed the output to 
                   a tar of text files. Now, the script does not select matches 
                   within a sensor-specific distance (which can now be done 
                   downstream), but grabs all matches within a maximum limit 
                   (nominally 30 km).
0.2.0    161102    Transitioned to a more generalized comparison. Modified 
                   output to be independent of stream. Also removed criterion 
                   of number of gauges and swath limitation in sounders; these 
                   can be implemented downstream.
0.1.1    160809    Added settings for ATMS.
0.1.0    160630    Code now can match a gauge to multiple GPROF pixels from the
                   same overpass. Added a selection of swath (for MHS). Removed 
                   use of masked array which is redundant anyway (masked value
                   of -9999 eliminates possibility of minimum distance).
0.0.5    160627    Added settings for MHS.
0.0.4    160520    Added settings for SSMIS. Changed the standard GPROF name 
                   from "GPROF1.4" to "V03". Updated paths.
0.0.3    160510    Further generalized for GPROF as well. Also switched to 
                   new GAG format.
0.0.2    160506    Generalized the script to use for other sensors.
0.0.1    160421    Created this script from iPython notebook.
'''
