#!/usr/bin/python

'''
Extracting satellite and gauge values from matched footprints.
Version 0.3.0

This script takes the matched footprints (from match_footprint-gauge.py) and 
uses it to extract precipitation-related data from GPROF/DPR and the gauges. 
The output contains text files of these values.

Jackson Tan (jackson.tan@nasa.gov)
2017 08 11
'''

import numpy as np
import h5py
import tarfile
import gzip
import sys
import os
from glob import glob
from calendar import monthrange
from datetime import datetime, timedelta


#--- Preliminaries ---#


year0, month0, year1, month1 = np.int32(sys.argv[1 : 5])
sensor, network = sys.argv[5 : 7]
try:
    dataVersion = sys.argv[7]
except IndexError:
    dataVersion = 'V05'
try:
    lag = int(sys.argv[8])
except IndexError:
    lag = 0

radars = ('Ku', 'DPR', 'DPRGMI')

datapath = '/home/jackson/Work/Data/footprint-gauge/'
gaugepath = '/media/jackson/Vault/'
if sensor in radars:
    satpath = '/media/jackson/Vault/DPR/%s/' % dataVersion
else:
    satpath = '/media/jackson/GPROF/%s/' % dataVersion

scriptVersion = '0.3.0'

nday = (datetime(year1, month1, monthrange(year1, month1)[1]) - 
        datetime(year0, month0, 1)).days + 1
platforms = {'Ku': ['GPM'], 'DPR': ['GPM'], 'DPRGMI': ['GPM'], 
             'GMI': ['GPM'], 'TMI': ['TRMM'], 'ATMS': ['NPP'], 
             'AMSR2': ['GCOMW1'], 'SSMIS': ['F16', 'F17', 'F18', 'F19'], 
             'MHS': ['METOPA', 'METOPB', 'NOAA18', 'NOAA19']}
sensors = platforms.keys()


#--- Settings ---#


mingauge = {s: {'PCMK': 6, 'WEGN': 12, 'OLPX': 3, 'WGEW': 12} 
            if s not in radars else
            {'PCMK': 3, 'WEGN': 6, 'OLPX': 2, 'WGEW': 6} for s in sensors}
maxdist = {s: 7.5 if s not in radars else 2.5 for s in sensors}


#--- Defining the required functions ---#


def getMonths(year, year0, month0, year1, month1):
    if   year0 == year1: return range(month0, month1 + 1)
    elif year  == year0: return range(month0, 13)
    elif year  == year1: return range(1, month1 + 1)
    else               : return range(1, 13)

def unpackMatchups(path, filestamp):
    files = ('footprints', 'gauges', 'distances')
    with tarfile.open('%s%s.tgz' % (path, filestamp)) as f:
        data = [f.extractfile('%s_%s.txt' % (filestamp, file)).read()\
                .decode().split('\n') for file in files]
    return data

def getFilename(date, platform, sensor, orbit):
    files = []
    for filename in glob('%s%s/%s/%s/*' % (satpath, date[0], date[1], date[2])):
        if ((('%s.%s-' % (platform, sensor) in filename) or
             ('%s.%s.' % (platform, sensor) in filename)) and
            (orbit in filename)):
            files.append(filename)
    if len(files):
        return files[0]
    else:
        return None


#--- Read and sort the matches ---#


mF = {platform: [] for platform in platforms[sensor]}
mG = {platform: [] for platform in platforms[sensor]}
mD = {platform: [] for platform in platforms[sensor]}

for platform in platforms[sensor]:
    for year in range(year0, year1 + 1):
        for month in getMonths(year, year0, month0, year1, month1):
            filestamp = '%s.%s_%4d%02d' % (platform, sensor, year, month)
            f, g, d = unpackMatchups('%s%s/' % (datapath, network), filestamp)
            mF[platform] += f
            mG[platform] += g
            mD[platform] += d


#--- Filter the matches for minimum gauges in a maximum distance ---#


fG = {platform: [] for platform in platforms[sensor]}
fD = {platform: [] for platform in platforms[sensor]}
fO = {platform: [] for platform in platforms[sensor]}
fI = {platform: [] for platform in platforms[sensor]}

for platform in platforms[sensor]:
    for F, G, D in zip(mF[platform], mG[platform], mD[platform]):
        if (sum([float(d) <= maxdist[sensor] for d in D.split()]) 
            >= mingauge[sensor][network]):
            fil = [dd for dd, d in enumerate(D.split()) 
                   if float(d) <= maxdist[sensor]]
            fG[platform].append([int(G.split()[ff]) for ff in fil])
            fD[platform].append([ff for ff in F.split()[:3]])
            fO[platform].append(F.split()[3])
            fI[platform].append((int(F.split()[4]), int(F.split()[5])))


#--- Read the satellite data ---#


P_sat = {platform: [] for platform in platforms[sensor]}    # precip rate
T_sat = {platform: [] for platform in platforms[sensor]}    # overpass time
Q_sat = {platform: [] for platform in platforms[sensor]}    # quality flag
S_sat = {platform: [] for platform in platforms[sensor]}    # surface index
A1_sat = {platform: [] for platform in platforms[sensor]}   # ancillary data
A2_sat = {platform: [] for platform in platforms[sensor]}   # ancillary data
I_sat = {platform: [] for platform in platforms[sensor]}   # indices
missingFiles = []

for platform in platforms[sensor]:
    for D, O, I in zip(fD[platform], fO[platform], fI[platform]):
        gproffile = getFilename(D, platform, sensor, O)
        if gproffile == None:
            missingFiles.append('%s.%s %s%s%s %s' % (platform, sensor, D[0], 
                                                     D[1], D[2], O))
            continue
        with h5py.File(gproffile, 'r') as f:

            if sensor in ('Ku', 'DPR'):
                P_sat[platform].append(f['NS/SLV/precipRateESurface'][I])
                Q_sat[platform].append(f['NS/FLG/qualityFlag'][I])
                S_sat[platform].append(f['NS/PRE/landSurfaceType'][I])
                A1_sat[platform].append(f['NS/SLV/precipRateNearSurface'][I])
                typePrecip = max(f['NS/CSF/typePrecip'][I], 0) // 10000000
                A2_sat[platform].append(typePrecip)
                basedir = 'NS'
            elif sensor == 'DPRGMI':
                P_sat[platform].append(f['NS/surfPrecipTotRate'][I])
                Q_sat[platform].append(-99)    # no quality flag for DPRGMI
                S_sat[platform].append(f['NS/Input/surfaceType'][I])
                A1_sat[platform].append(f['NS/surfPrecipTotRateSigma'][I])
                typePrecip = (max(f['NS/Input/precipitationType'][I], 0) // 
                              10000000)
                A2_sat[platform].append(typePrecip)
                basedir = 'NS'
            else:
                P_sat[platform].append(f['S1/surfacePrecipitation'][I])
                Q_sat[platform].append(f['S1/qualityFlag'][I])
                S_sat[platform].append(f['S1/surfaceTypeIndex'][I])
                A1_sat[platform].append(f['S1/probabilityOfPrecip'][I])
                if dataVersion == 'V04':
                    A2_sat[platform].append(f['S1/convectPrecipFraction'][I])
                else:
                    A2_sat[platform].append(f['S1/convectivePrecipitation'][I])
                basedir = 'S1'

            year = f['%s/ScanTime/Year' % basedir][I[0]]
            month = f['%s/ScanTime/Month' % basedir][I[0]]
            day = f['%s/ScanTime/DayOfMonth' % basedir][I[0]]
            hour = f['%s/ScanTime/Hour' % basedir][I[0]]
            minute = f['%s/ScanTime/Minute' % basedir][I[0]]
            second = f['%s/ScanTime/Second' % basedir][I[0]]
            OPtime = (datetime(year, month, day, hour, minute, second) + 
                      timedelta(minutes = lag))
            T_sat[platform].append(OPtime)
        I_sat[platform].append(I)


#--- Read the gauge data ---#


if network == 'PCMK':

    pcmkpath = '%sGAG/PCMK/' % gaugepath

    # get snow days
    snowdays = []
    with open('%ssnowdays' % pcmkpath, 'r') as f:
        for day in f.readlines():
            year, month, day = day.split()
            snowdays.append(datetime(int(year), int(month), int(day)))
            snowdays.append((datetime(int(year), int(month), int(day)) + 
                             timedelta(days = 1)))
    snowdays = list(set(snowdays))

    # add missing data in PCMK to the "snowdays" list for masking
    days = (datetime(2015, 10, 7) - datetime(2015, 9, 21)).days + 1
    missing = [(datetime(2015, 9, 21) + timedelta(days = day))
               for day in range(days)]
    snowdays += missing

    # read days for PCMK values that have been flagged for removal
    flagged = {}
    with tarfile.open('%sfilters.tgz' % (pcmkpath,), 'r:gz') as f:
        for gauge in range(1, 26):
            d1 = str(f.extractfile('moved_NASA%04d.txt' % gauge).read(), 
                     'utf-8').split('\n')[3 : -1]
            d2 = str(f.extractfile('pair_NASA%04d.txt' % gauge).read(), 
                     'utf-8').split('\n')[3 : -1]
            d3 = str(f.extractfile('next_NASA%04d.txt' % gauge).read(), 
                     'utf-8').split('\n')[3 : -1]
            filter2 = [datetime(int(d[0 : 4]), int(d[6 : 8]), int(d[10 : 12])) 
                       for d in d2]
            filter3 = [datetime(int(d[0 : 4]), int(d[6 : 8]), int(d[10 : 12])) 
                       for d in d3]
            filter1 = []
            for ii in range(len(d1) // 3):
                m = d1[ii * 3 + 1]
                date1 = datetime(int(m[0 : 4]), int(m[6 : 8]), int(m[10 : 12]))
                m = d1[ii * 3 + 2]
                date2 = datetime(int(m[0 : 4]), int(m[6 : 8]), int(m[10 : 12]))
                filter1 += [date1 + timedelta(days = n) 
                            for n in range((date2 - date1).days + 1)]
            flagged[gauge] = sorted(list(set(filter1 + filter2 + filter3)))

    # read the raw PCMK data
    data = {}
    for gg, gauge in enumerate(range(1, 26)):
        for bb, bucket in enumerate(('A', 'B')):
            file = '%scombined/IOWA-NASA%04d_%s.gag' % (pcmkpath, gauge, bucket)
            with open(file, 'r') as f:
                data[gauge, bucket] = [l.split() for l in f.readlines()[2:]]

    # convert to 1-min accumulations
    P_pcmk = {}
    t0 = datetime(year0, month0, 1)
    t1 = datetime(year1, month1, monthrange(year1, month1)[1]) + \
                  timedelta(days = 1)
    for gauge in range(1, 26):
        P_pcmk[gauge] = np.zeros([2, nday * 24 * 60], dtype = np.float32)
        for bb, bucket in enumerate(('A', 'B')):
            for ll in data[gauge, bucket]:
                tiptime = datetime(int(ll[0]), int(ll[1]), int(ll[2]), 
                                   int(ll[4]), int(ll[5]), int(ll[6]))
                if tiptime >= t0 and tiptime < t1:
                    minute = int((tiptime - t0).total_seconds()) // 60
                    P_pcmk[gauge][bb, minute] += float(ll[7]) * 60

    # mask 1-min accumulation for snow days and flagged days
    mask = np.zeros([2, nday * 24 * 60], dtype = np.bool)
    for snowday in snowdays:
        if t0 <= snowday < t1:
            index = int((snowday - t0).total_seconds()) // 60
            mask[:, index : index + 24 * 60 * 2] = True
    for gauge in range(1, 26):
        P_pcmk[gauge] = np.ma.masked_where(mask, P_pcmk[gauge])
    for gauge in range(1, 26):
        mask = np.zeros([2, nday * 24 * 60], dtype = np.bool)
        for flaggedday in flagged[gauge]:
            if t0 <= flaggedday < t1:
                index = int((flaggedday - t0).total_seconds()) // 60
                mask[:, index : index + 24 * 60] = True
        P_pcmk[gauge] = np.ma.masked_where(mask, P_pcmk[gauge])

    # extract gauge values over OP time
    nmin = 2    # number of minutes in both directions to extract
    P_gauge = {}
    for platform in platforms[sensor]:
        nmatch = len(fO[platform]) - len(missingFiles)
        P_gauge[platform] = [[] for _ in range(nmatch)]
        for mm in range(nmatch):
            if T_sat[platform][mm] + timedelta(seconds = nmin * 60) >= t1:
                continue    # skip if OP time extends beyond period
            t = int((T_sat[platform][mm] - t0).total_seconds()) // 60
            P = np.ma.array([P_pcmk[gg][:, t - nmin : t + nmin + 1] 
                             for gg in fG[platform][mm]])
            if np.ma.count(P):
                P_gauge[platform][mm] += np.mean(P, 2).compressed().tolist()

if network == 'WEGN':

    wegnpath = '%sWegenerNet/Station/' % gaugepath
    datestamp = '2014-03-31d00h00m_2017-04-01d23h55m'

    # read the raw WEGN data
    filename = 'WN_L2_V6_BD_St1_%s_UTC.csv.gz' % datestamp
    T_wegn = []
    with gzip.open('%s%s' % (wegnpath, filename), 'r') as f:
        for rawline in f.readlines()[1:]:
            line = rawline.decode().split(',')
            T_wegn.append(datetime(int(line[1][:4]), int(line[1][5 : 7]), 
                                   int(line[1][8 : 10]), int(line[1][11 : 13]), 
                                   int(line[1][14 : 16]), int(line[1][17 : 19]))
                          - timedelta(minutes = 5))
    T_wegn = np.array(T_wegn)

    # truncate to desired time period
    t0 = datetime(year0, month0, 1)
    t1 = datetime(year1, month1, monthrange(year1, month1)[1]) + \
                  timedelta(days = 1)
    i0 = np.where(T_wegn == t0)[0][0]
    i1 = np.where(T_wegn == t1)[0][0]
    T_wegn = T_wegn[i0 : i1]

    # read the 5-min accumulations
    P_wegn = {}
    for stn in range(1, 152):
        filename = 'WN_L2_V6_BD_St%d_%s_UTC.csv.gz' % (stn, datestamp)
        value, qflag = [], []
        with gzip.open('%s%s' % (wegnpath, filename), 'r') as f:
            for rawline in f.readlines()[1:]:
                line = rawline.decode().split(',')
                if int(line[0]) != stn:
                    sys.exit('Error: station number does not match CSV record.')
                value.append(float(line[2]))
                qflag.append(int(line[3]))
        P_wegn[stn] = (np.ma.masked_where(qflag[i0 : i1], value[i0 : i1]) * 12)

    # extract gauge values over OP time
    nts = 0    # number of WEG timesteps (5 min) to extract in both directions
    P_gauge = {}
    for platform in platforms[sensor]:
        nmatch = len(fO[platform]) - len(missingFiles)
        P_gauge[platform] = [[] for _ in range(nmatch)]
        for mm in range(nmatch):
            if T_sat[platform][mm] + timedelta(minutes = nts * 5) >= t1:
                continue    # skip if OP time extends beyond period
            OPtime = T_sat[platform][mm] - t0
            t = int(OPtime.total_seconds() / 60 / 5)
            P = np.ma.array([P_wegn[gg][t - nts : t + nts + 1] 
                             for gg in fG[platform][mm]])
            if np.ma.count(P):
                P_gauge[platform][mm] += [ii for ii in np.mean(P, 1).tolist()
                                          if ii != None]

if network == 'OLPX':

    olpxpath = '%sGAG/OLPX/combined/' % gaugepath

    gauges = [26, 27, 28, 29, 31, 32, 33, 34, 35, 37, 39, 40, 
              41, 42, 43, 44, 45, 101]

    # read the raw OLPX data
    data = {}
    for gg, gauge in enumerate(gauges):
        for bb, bucket in enumerate(('A', 'B')):
            file = '%sOLYMPEx-NASA%04d_%s.gag' % (olpxpath, gauge, bucket)
            with open(file, 'r') as f:
                data[gauge, bucket] = [l.split() for l in f.readlines()[2:]]

    # convert to 1-min accumulations
    P_olpx = {}
    t0 = datetime(year0, month0, 1)
    t1 = datetime(year1, month1, monthrange(year1, month1)[1]) + \
                  timedelta(days = 1)
    for gauge in gauges:
        P_olpx[gauge] = np.zeros([2, nday * 24 * 60], dtype = np.float32)
        for bb, bucket in enumerate(('A', 'B')):
            for ll in data[gauge, bucket]:
                tiptime = datetime(int(ll[0]), int(ll[1]), int(ll[2]), 
                                   int(ll[4]), int(ll[5]), int(ll[6]))
                if tiptime >= t0 and tiptime < t1:
                    minute = int((tiptime - t0).total_seconds()) // 60
                    P_olpx[gauge][bb, minute] += float(ll[7]) * 60

    # extract gauge values over OP time
    nmin = 2    # number of minutes in both directions to extract
    P_gauge = {}
    for platform in platforms[sensor]:
        nmatch = len(fO[platform]) - len(missingFiles)
        P_gauge[platform] = [[] for _ in range(nmatch)]
        for mm in range(nmatch):
            if T_sat[platform][mm] + timedelta(seconds = nmin * 60) >= t1:
                continue    # skip if OP time extends beyond period
            t = int((T_sat[platform][mm] - t0).total_seconds()) // 60
            P = np.ma.array([P_olpx[gg][:, t - nmin : t + nmin + 1] 
                             for gg in fG[platform][mm]])
            if np.ma.count(P):
                P_gauge[platform][mm] += np.mean(P, 2).flatten().tolist()

if network == 'WGEW':

    wgewpath = '%sWGEW/' % gaugepath

    td2min = lambda x : int(x.total_seconds() / 60)

    # initialize the number of minutes in the chosen period
    t0 = datetime(year0, month0, 1)
    t1 = datetime(year1, month1, monthrange(year1, month1)[1]) + \
                  timedelta(days = 1)
    minInPeriod = int((t1 - t0).total_seconds() / 60)

    # read the raw WGEW data into 1-min accumulations
    P_wgew = {}
    for file in sorted(glob('%sRainGaugeData/wgew_*.txt' % wgewpath)):

        with open(file, 'r') as f:
            data = [ll.split(',') for ll in f.readlines()[3:]]
        gauge = int(data[0][0])
        P_wgew[gauge] = np.zeros(minInPeriod, dtype = 'f4')

        for line in data:
            if float(line[3]) == 0:
                tevent = (datetime.strptime(line[1] + line[2], '%m/%d/%Y%H:%M')
                          + timedelta(hours = 7))    # convert from MST to UTC
                telapsed = timedelta(0)
                p = float(line[5])
            else:
                ts = tevent + telapsed
                telapsed = timedelta(minutes = float(line[3]))
                te = tevent + telapsed
                if   te > t0 and ts < t1:
                    i0 = max(td2min(ts - t0), 0)
                    i1 = min(td2min(te - t0), minInPeriod)
                    P_wgew[gauge][i0 : i1] = p
                p = float(line[5])

    # mask gauges 101-109 prior to start of operation on 2015 07 07
    t = int((datetime(2015, 7, 8) - t0).total_seconds() / 60)
    if t > 0:
        for gauge in range(101, 110):
            P_wgew[gauge][:t] = np.nan

    ## mask days below freezing (in development)
    #with open('%sMeteorologicalSiteData/km14.out' % wgewpath, 'r') as f:
    #    data = [ii.split() for ii in f.readlines()]
    #Tday = [[] for _ in range(365)]
    #for line in data:
    #    Tday[int(line[3]) - 1].append(float(line[5]))
    #belowFreezing = [any([ii < 0 for ii in tt]) for tt in Tday]

    # extract gauge values over OP time
    nmin = 2    # number of minutes in both directions to extract
    P_gauge = {}
    for platform in platforms[sensor]:
        nmatch = len(fO[platform])
        P_gauge[platform] = [[] for _ in range(nmatch)]
        for mm in range(nmatch):
            if T_sat[platform][mm] + timedelta(seconds = nmin * 60) >= t1:
                continue    # skip if OP time extends beyond period
            t = int((T_sat[platform][mm] - t0).total_seconds()) // 60
            P = np.array([P_wgew[gg][t - nmin : t + nmin + 1] 
                          for gg in fG[platform][mm]])
            P_gauge[platform][mm] += [ii for ii in np.mean(P, 1).tolist()
                                      if not np.isnan(ii)]


#--- Filter out masked values and collapse the platform dimension ---#


X1, X2, Y, Q, S, T, A1, A2, P, I = [], [], [], [], [], [], [], [], [], []

for platform in platforms[sensor]:
    nmatch = len(fO[platform]) - len(missingFiles)
    mask1 = [P_sat[platform][mm] < 0 for mm in range(nmatch)]
    mask2 = [len(ii) == 0 for ii in P_gauge[platform]]
    mask3 = [len(ii) < mingauge[sensor][network] for ii in P_gauge[platform]]
    mask = np.array(mask1) + np.array(mask2) + np.array(mask3)
    X1 += [np.mean(P_gauge[platform][mm], dtype = np.float32) 
           for mm in range(nmatch) if not mask[mm]]
    X2 += [P_gauge[platform][mm] for mm in range(nmatch) if not mask[mm]]
    Y  += [P_sat[platform][mm] for mm in range(nmatch) if not mask[mm]]
    Q  += [Q_sat[platform][mm] for mm in range(nmatch) if not mask[mm]]
    S  += [S_sat[platform][mm] for mm in range(nmatch) if not mask[mm]]
    A1 += [A1_sat[platform][mm] for mm in range(nmatch) if not mask[mm]]
    A2 += [A2_sat[platform][mm] for mm in range(nmatch) if not mask[mm]]
    T  += [(T_sat[platform][mm] - timedelta(minutes = lag))
           .strftime('%y%m%d %H%M%S') for mm in range(nmatch) if not mask[mm]]
    P  += [platform for mm in range(nmatch) if not mask[mm]]
    I  += [I_sat[platform][mm] for mm in range(nmatch) if not mask[mm]]


#--- Write the values to file ---#


# set the directory
outpath = '%sMatchedValues/%s/' % (datapath, dataVersion)
if not os.path.exists(outpath): os.makedirs(outpath)
daterange = '%4d%02d-%4d%02d' % (year0, month0, year1, month1)

# filenames corresponding to the variables
file01 = {s: 'gaugePrecipitation.txt' for s in sensors}
file02 = {s: 'gaugePrecipitationAll.txt' for s in sensors}
file03 = {s: 'surfacePrecipitation.txt' if s not in radars
          else 'surfPrecipTotRate.txt' if s == 'DPRGMI'
          else 'precipRateESurface.txt' for s in sensors}
file04 = {s: 'overpassTime.txt' for s in sensors}
file05 = {s: 'probabilityOfPrecip.txt' if s not in radars
          else 'surfPrecipTotRateSigma.txt' if s == 'DPRGMI'
          else 'precipRateNearSurface.txt' for s in sensors}
file06 = {s: 'convectivePrecipitation.txt' if s not in radars
          else 'typePrecip.txt' for s in sensors}
file07 = {s: 'surfaceTypeIndex.txt' if s not in radars
          else 'surfaceType.txt' if s == 'DPRGMI'
          else 'landSurfaceType.txt' for s in sensors}
file08 = {s: 'qualityFlag.txt' for s in sensors}
file09 = {s: 'indices.txt' for s in sensors}
file10 = {s: 'platform.txt' for s in sensors}
file11 = {s: 'version.txt' for s in sensors}
file12 = {s: 'missingFiles.txt' for s in sensors}

# append random string to temp files to avoid overwriting from parallel procs
tmp = np.random.random()

# write the text files
with open('%s%s%s' % (outpath, file01[sensor], tmp), 'w') as f:
    for ii in X1:
        f.write('%s\n' % ii)
with open('%s%s%s' % (outpath, file02[sensor], tmp), 'w') as f:
    for ii in X2:
        fmt = '{:7.3f} ' * len(ii) + '\n'
        f.write(fmt.format(*ii))
with open('%s%s%s' % (outpath, file03[sensor], tmp), 'w') as f:
    for ii in Y:
        f.write('%s\n' % ii)
with open('%s%s%s' % (outpath, file04[sensor], tmp), 'w') as f:
    for ii in T:
        f.write('%s\n' % ii)
with open('%s%s%s' % (outpath, file05[sensor], tmp), 'w') as f:
    for ii in A1:
        f.write('%s\n' % ii)
with open('%s%s%s' % (outpath, file06[sensor], tmp), 'w') as f:
    for ii in A2:
        f.write('%s\n' % ii)
with open('%s%s%s' % (outpath, file07[sensor], tmp), 'w') as f:
    for ii in S:
        f.write('%s\n' % ii)
with open('%s%s%s' % (outpath, file08[sensor], tmp), 'w') as f:
    for ii in Q:
        f.write('%s\n' % ii)
with open('%s%s%s' % (outpath, file09[sensor], tmp), 'w') as f:
    for ii in I:
        f.write('%4s  %3s\n' % (ii[0], ii[1]))
with open('%s%s%s' % (outpath, file10[sensor], tmp), 'w') as f:
    for ii in P:
        f.write('%s\n' % ii)
with open('%s%s%s' % (outpath, file11[sensor], tmp), 'w') as f:
    f.write('Produced by getvalues_footprint-gauge.py ver. %s' % scriptVersion)
    f.write(' at %s.' % datetime.utcnow().strftime('%Y/%m/%d %H:%M:%S UTC'))
with open('%s%s%s' % (outpath, file12[sensor], tmp), 'w') as f:
    for ii in set(missingFiles):
        f.write('%s\n' % ii)

# place the text files in a zipped tarball
filename = '%s-%s%+d_%s.tgz' % (sensor, network, lag, daterange)
with tarfile.open('%s%s' % (outpath, filename), 'w:gz') as f:
    for file in (file01, file02, file03, file04, file05, file06, file07, 
                 file08, file09, file10, file11, file12):
        f.add('%s%s%s' % (outpath, file[sensor], tmp), arcname = file[sensor])
        os.remove('%s%s%s' % (outpath, file[sensor], tmp))


'''
Version history
0.3.0    170811    Added more outputs: quality flag, surface index and indices.
0.2.10   170809    Fixed an error on the platform output.
0.2.9    170802    Masked WGEW gauges before their start of operation.
0.2.8    170720    Added QC filters for PCMK.
0.2.7    170707    Reduced accumulation period to 5-min for all networks; this
                   is so as to better reflect the instantaneous nature of the 
                   satellite estimates. The results from lag-accumulation 
                   should not be a factor.
0.2.6    170703    Added platform as another output.
0.2.5    170626    Fixed handling of mask values. Is now consistent with 
                   accum-lag_footprint-gauge.py.
0.2.4    170623    Corrected WGEW for MST to UTC. Also, trivial code cleanup.
0.2.3    170622    Fixed an error reading WGEW rain rates.
0.2.2    170621    Added WGEW network.
0.2.1    170602    Improved handling of GPROF files, especially for missing 
                   files. Also increased accumulation period to HSS maximum.
0.2.0    170525    Modified code to work with new WEGN data format.
0.1.16   170417    Modified input variables.
0.1.15   170330    Fixed a minor error for DPRGMI due to duplicated variable;
                   changed it to an alternative variable.
0.1.14   170329    Added platforms for GPM.DPR and GPM.DPRGMI. Also generalized
                   Ku to multiple options.
0.1.13   170302    Fixed an error in handling snowday/missing mask in PCMK.
0.1.12   170301    Fixed stupid typo "gaugePreciptiationAll".
0.1.11   170227    Added another output 'gaugePrecipitationAll'. Also fixed an 
                   error involving removing masked elements (affects WEGN).
0.1.10   170223    Reverted the previous accumulation settings based on new PSS
                   analysis.
0.1.9    170222    Modified accumulation period to the average of accumulation 
                   above the 95th percentile of PSS scores as a function of
                   threshold-lag-accum. Also improved handling of variable 
                   names between V04 and V05.
0.1.8    170103    Updated variables names to those used in GPROF V05.
0.1.7    161228    Added control over GPROF/DPR version (both in input and in
                   output).
0.1.6    161219    Fixed an error in overpassTime when lag is nonzero. Added 
                   version output. Also reduced mingauge for OLPX.
0.1.5    161208    Added convectPrecipFraction as another ancillary output for
                   the GPROF estimates. Also improved output.
0.1.4    161202    Fixed an error in maxdist; should be 2.5 km for Ku (no 
                   effect for GPROF). Also adjusted mingauge correspondingly.
0.1.3    161128    Doubled the minimum number of gauges based on feedback.
0.1.2    161125    Fixed an error in input filename in WEGN. Implemented WEGN
                   correction factors. Added OP time into output. 
0.1.1    161123    Added safeguard against OP time exceeding record time. Added
                   setting for lag between OP time and gauges. Reduced PCMK and 
                   OLPX gauges to minute resolution for memory efficiency. Also 
                   fixed minor errors in accumulation period for WEGN.
0.1.0    161122    Integrated DPR Ku (match_Ku-gauge.py). Added OLPX.
0.0.2    161121    Moved masking of snowdays for PCMK to the gauges.
0.0.1    161110    Created this script from iPython notebook.
'''
