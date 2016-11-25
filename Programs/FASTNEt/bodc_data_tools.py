# -*- coding: utf-8 -*-
"""
bodc_data_tools.py

Collection of functions to process BODC data. Initially at least the aim is to read in FASTNEt mooring data for processing.

Created on Fri 4 Dec 2015

@author: jeff
"""
import numpy as np
#import matplotlib.pyplot as plt
#import scipy.optimize as optimization
#import datetime
#import re
#import bisect
#import requests,json
import csv
from itertools import islice
from scipy.io import netcdf
import numpy.ma as ma

#==============================================================================
#==============================================================================
def readBodcDataGivenLoc(mooring, dirpath):
    """
    Function to return the temperature, time and depth along a mooring
    Each file correspond to an individual instrument bodc NetCDF files: e.g. b1198294.qxf
    The lat and lon of each mooring is passed to this routine.
    The routine finds all the instruments at this location, and extracts, temperature and time
    These temperature series are then ordered by depth
    A mask is applied to the outlier/missing data
    Temperature, depth and time are returned
    Input:
        lonlatpair : location of the mooring. e.g. [('-9.14678', '48.62975')]
        path where to find the bocd netcdf files
    Output:
        temperature (2D array)
        depth (1D array)
        time (2D array)
    Usage:
        [temperature,depth,time] = readBodcDataGivenLoc(mooring, CRUISE_DIR)
        temperature = {}
        time_depth={}
        max = 0
    """
    temperature = {}
    time_depth={}
    max = 0
    for instru in mooring:

        #load temperature for each depth
        filepath = dirpath + 'b' + instru['BODC reference'] + '.qxf'
        f = netcdf.netcdf_file(filepath, 'r')
        temp = f.variables['TEMPPR01'].data
        depth = instru['Series depth maximum (m)']
        if temperature.has_key(depth):
            temperature[depth].append(temp)
        else:
            temperature[depth] = [temp]
        if max < len(temp):
            max = len(temp)

        #load time for each depth
        day = f.variables['AADYAA01'].data
        hour = f.variables['AAFDZZ01'].data
        time = (day+hour)
        if time_depth.has_key(depth):
            time_depth[depth].append(time)
        else:
            time_depth[depth] = [time]
        f.close()

    # order outputs by increasing depth
    nbr_depths = len(temperature.keys())
    ordered_temperature = np.zeros((nbr_depths, max))
    ordered_time = np.zeros((nbr_depths, max))
    ordered_depth = np.zeros((nbr_depths))
    for index, depth in enumerate(sorted(temperature, key=float)):
        n = np.array(temperature[depth])
        ordered_temperature[index,0:n.shape[1]] = n
        n = np.array(time_depth[depth])
        ordered_time[index,0:n.shape[1]] = n
        ordered_depth[index] = depth


    # Mask outliers (0 and val outside (mean +- 1.8std))
    ordered_temperature_masked = ma.masked_where(ordered_temperature == 0,ordered_temperature)
    ordered_time_masked = ma.masked_where(ordered_time == 0,ordered_time)
    for index in range(len(ordered_depth)):
        ampl_valid = 1.8
        temp_std = np.std(ordered_temperature_masked)
        temp_min = ordered_temperature_masked[index].mean() - ampl_valid * temp_std
        temp_max = ordered_temperature_masked[index].mean() + ampl_valid * temp_std
        ordered_temperature_masked[index] = ma.masked_outside(ordered_temperature_masked[index],\
                                                              temp_min, \
                                                              temp_max)
    ordered_time_masked = ma.masked_where(ordered_temperature_masked == False, ordered_time_masked)

    return [ordered_temperature_masked, ordered_time_masked, ordered_depth]

def readBodcCSV(csv_filepath):
    """
    Reads in metadata from a bodc summary file
    These data are needed since the actual data files do not contain it.
    In particular this function will extract the depth, lat and lon of the instrument.
    NB each BODC file contains a timeseries for a single instrument. Therefore profiles need to be reconstructed from lat,lon pairings

    Input:
        None
    Output:
        resultdictionary
        result[0]['BODC reference'] : instrument reference identifiers - e.g. 1198859
        result[0].keys() : displays available keys
        CalledFromAbove : now set to CalledFromAbove=1
    Usage:
        resultdictionary = readBodcData()
        print resultdictionary[0]['BODC reference']
        print resultdictionary[0].keys()
    """
    with open(csv_filepath, "r") as f:
        '''Pass preamble to find where the data starts'''
        n = 0
        for line in f.readlines():
            n += 1
            if 'BODC reference,Oceanographic data type' in line: # line with field names was found
                h = line.split(',')
                break
        f.close()
    # read in the data from the desired line with the appropriate header labels
    f = islice(open(csv_filepath, "r"), n, None)

    content = []
    reader = csv.DictReader(f, fieldnames = h)
    for row in reader:
        content.append(row);
    return content


def readBodcInterpTime(temperature,time,freq):
    """
    Interpolates the temperature along a mooring on a regular time axis
     It is needed as each instrument along the mooring has its own sampling frequency.

      Input:

      Output:
      Usage:
    """
    time_start = int(time.min())
    time_end = int(time.max())+1

    tot_pts=(time_end - time_start) * freq
    print('Interpolation on {} days, i.e. {} points'.format(time_end-time_start, tot_pts))
    freq_taxis = np.linspace(time_start, time_end, num=tot_pts)
    freq_temperature = np.zeros((len(temperature),len(freq_taxis)))

    loop_max = len(temperature)
    for index in range(loop_max):
        time_ind = time[index]
        temperature_ind = temperature[index]
        freq_temperature[index] = np.interp(freq_taxis,\
                                              time[index].compressed(),\
                                              temperature[index].compressed(),\
                                              left=0, right=0)
    freq_temperature = ma.masked_where(freq_temperature == 0, freq_temperature)
    return[freq_temperature,freq_taxis]


def readBodcWriteNC(ncname,temperature,time_array,depth_array,longitude,latitude):
    """
    Writes a netcdf file with temperature function of depth and time
    """
    f = netcdf.netcdf_file(ncname, 'w')
    f.history = 'Temperature along FASTNEt mooring, computed from bodc files'

    f.createDimension('time',len(time_array))
    f.createDimension('depth',len(depth_array))
    f.createDimension('coord',1)

    time = f.createVariable('time', 'f', ('time',))
    time[:] = time_array
    time.units = 'days since 1760-01-01'

    depth = f.createVariable('depth','f',('depth',))
    depth[:] = depth_array

    lon = f.createVariable('lon','f',('coord',))
    lon[0] = float(longitude)
    lat = f.createVariable('lat','f',('coord',))
    lat[0] = float(latitude)

    tempobs = f.createVariable('Temp_Obs','d',('depth','time',))
    tempobs[:] = temperature

    f.close()


#==============================================================================
def readBodcMooringTemp(mooring, filedir):
    """
    Function to return the temperature, time and depth along a mooring
    Each file correspond to an individual instrument bodc NetCDF files: e.g. b1198294.qxf
    The lat and lon of each mooring is passed to this routine.
    The routine finds all the instruments at this location, and extracts, temperature and time
    These temperature series are then ordered by depth
    A mask is applied to the outlier/missing data
    Temperature, depth and time are returned
    Input:
        lonlatpair : location of the mooring. e.g. [('-9.14678', '48.62975')]
        path where to find the bocd netcdf files
    Output:
        temperature (2D array)
        depth (1D array)
        time (2D array)
    Usage:
        [temperature,depth,time] = readBodcMooringTemp(mooring, CRUISE_DIR)
        temperature = {}
        time_depth={}
        max = 0
    """
    filepath_obs = filedir + 'mooring_' + mooring + '.nc'
    f = netcdf.netcdf_file(filepath_obs, 'r')
    temp_obs = f.variables['Temp_Obs'].data
    time_obs = f.variables['time'].data
    depth_obs = f.variables['depth'].data
    f.close()

    dum = temp_obs.copy()
    temp_obs = dum

    temp_obs[temp_obs==0] = np.NaN
    temp_obs = ma.array(temp_obs,mask=np.isnan(temp_obs))

    return [temp_obs,time_obs,depth_obs]

