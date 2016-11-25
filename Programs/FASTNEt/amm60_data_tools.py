# -*- coding: utf-8 -*-
"""
amm60_data_tools.py

Collection of functions to process AMM60 and obs data.
Created on Fri 4 Dec 2015

@author: jeff
"""
import numpy as np
#import matplotlib.pyplot as plt
#import scipy.optimize as optimization
import datetime
#import re
#import bisect
#import requests,json
import csv
from itertools import islice
from scipy.io import netcdf
import numpy.ma as ma
import h5py
from jdcal import gcal2jd, jd2gcal
#==============================================================================
#==============================================================================
def writeNCmodelobs(ncname,time_array,depth_array,longitude,latitude,tmod,tobs):
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
    time.units = 'days since 1950-01-01'

    depth = f.createVariable('depth','f',('depth',))
    depth[:] = depth_array

    lon = f.createVariable('lon','f',('coord',))
    lon[0] = float(longitude)
    lat = f.createVariable('lat','f',('coord',))
    lat[0] = float(latitude)

    temp_obs = f.createVariable('temp_obs','d',('depth','time',))
    temp_obs[:,:] = tobs

    temp_mod = f.createVariable('temp_mod','d',('depth','time',))
    temp_mod[:,:] = tmod

    f.close()


def readNCmooring_hdf5(filepath,st,en):
    """
    Read the NEMO outputs (mooring time series)
    """

    f = h5py.File(filepath, 'r')
    temp_mod = f['thetao'][:]
    time_mod = f['time_counter'][:]
    depth_mod = f['depth'][:]
    lon = f['nav_lon'][:]
    lat = f['nav_lat'][:]
    f.close()

    indtime_st=st
    indtime_en=en

    longitude = lon.mean()
    latitude = lat.mean()

    time_mod = time_mod[indtime_st:indtime_en]

    depth_mod = depth_mod.mean(axis=3)
    depth_mod = depth_mod.mean(axis=2)
    depth_mod = np.transpose(depth_mod)
    depth_mod = depth_mod[:,indtime_st:indtime_en]

    temp_mod = temp_mod.mean(axis=3)
    temp_mod = temp_mod.mean(axis=2)
    temp_mod = ma.masked_where(temp_mod == 0, temp_mod)
    temp_mod = np.transpose(temp_mod)
    temp_mod = temp_mod[:,indtime_st:indtime_en]

    dum = temp_mod.copy()
    temp_mod = dum

    temp_mod[temp_mod==0] = np.NaN
    temp_mod = ma.array(temp_mod,mask=np.isnan(temp_mod))

    return [temp_mod,time_mod,depth_mod,latitude,longitude]

def readNCmooring_nna(filepath,st,en):
    """
    Read the NEMO outputs (mooring time series)
    """

    f = h5py.File(filepath, 'r')
    temp_mod = f['votemper'][:]
    time_mod = f['time_counter'][:]
    depth_mod = f['deptht'][:]
    lon = f['nav_lon'][:]
    lat = f['nav_lat'][:]
    f.close()

    indtime_st=st
    indtime_en=en

    longitude = lon.mean()
    latitude = lat.mean()

    #time_mod = time_mod #[indtime_st:indtime_en]
    time_mod = time_mod[indtime_st:indtime_en]

    #depth_mod = depth_mod.mean(axis=3)
    #depth_mod = depth_mod.mean(axis=2)
    #depth_mod = np.transpose(depth_mod)
    #depth_mod = depth_mod[indtime_st:indtime_en]

    temp_mod = temp_mod.mean(axis=3)
    temp_mod = temp_mod.mean(axis=2)
    temp_mod = ma.masked_where(temp_mod == 0, temp_mod)
    temp_mod = np.transpose(temp_mod)
    temp_mod = temp_mod[:,indtime_st:indtime_en]

    dum = temp_mod.copy()
    temp_mod = dum

    temp_mod[temp_mod==0] = np.NaN
    temp_mod = ma.array(temp_mod,mask=np.isnan(temp_mod))

    return [temp_mod,time_mod,depth_mod,latitude,longitude]

def readNCmooring_nc(filepath,st,en):
    """
    Read the NEMO outputs (mooring time series)
    """

    f = netcdf.netcdf_file(filepath, 'r')
    depth_mod = f.variables['depth'].data
    temp_mod = f.variables['thetao'].data
    time_mod = f.variables['time_counter'].data
    lon = f.variables['nav_lon'].data
    lat = f.variables['nav_lat'].data
    f.close()

    indtime_st=st
    indtime_en=en

    longitude = lon.mean()
    latitude = lat.mean()

    time_mod = time_mod[indtime_st:indtime_en]

    depth_mod = depth_mod.mean(axis=3)
    depth_mod = depth_mod.mean(axis=2)
    depth_mod = np.transpose(depth_mod)
    depth_mod = depth_mod[:,indtime_st:indtime_en]

    temp_mod = temp_mod.mean(axis=3)
    temp_mod = temp_mod.mean(axis=2)
    temp_mod = ma.masked_where(temp_mod == 0, temp_mod)
    temp_mod = np.transpose(temp_mod)
    temp_mod = temp_mod[:,indtime_st:indtime_en]

    dum = temp_mod.copy()
    temp_mod = dum

    temp_mod[temp_mod==0] = np.NaN
    temp_mod = ma.array(temp_mod,mask=np.isnan(temp_mod))

    return [temp_mod,time_mod,depth_mod,longitude,latitude]


def interpMooringModelData(temp_mod,time_mod,depth_mod,temp_obs,time_obs,depth_obs):
    """
    # Interpolate obs data on model grid (depth, time)
    # 1. interpolate obs on model time (time_lev)
    # 2. interpolate mod and obs on mean model depth (depth_lev)
    """

    time_lev = time_mod/86400
    #time_delay = 365*190 + (190/4) 1
    #time_obs_1950 = time_obs - time_delay
    obs_st = gcal2jd(1760,1,1)[0]+gcal2jd(1760,1,1)[1]+0.5
    mod_st = gcal2jd(1950,1,1)[0]+gcal2jd(1950,1,1)[1]+0.5
    delta_cal = mod_st - obs_st
    time_obs_1950 = time_obs -delta_cal
    depth_lev = depth_mod.mean(axis=1)

    # 1.
    temp_obs_tinterp = np.zeros((len(depth_obs),len(time_lev)))
    for index in range(len(depth_obs)):
        temp_obs_tinterp[index,:] = np.interp(time_lev,\
                                          time_obs_1950,\
                                          temp_obs[index,:])

    # 2
    temp_obs_interp = np.zeros((len(depth_lev),len(time_lev)))
    temp_mod_interp = np.zeros((len(depth_lev),len(time_lev)))

    for index in range(len(time_lev)):
        temp_obs_interp[:,index] = np.interp(depth_lev,\
                                depth_obs,\
                                temp_obs_tinterp[:,index],\
                                left=None, right=None)

        temp_mod_interp[:,index] = np.interp(depth_lev,\
                                depth_mod[:,index],\
                                temp_mod[:,index])

    return [temp_mod_interp,temp_obs_interp,time_lev,depth_lev]

def interpMooringNNAData(temp_mod,time_mod,depth_mod,temp_obs,time_obs,depth_obs):
    """
    # Interpolate obs data on model grid (depth, time)
    # 1. interpolate obs on model time (time_lev)
    # 2. interpolate mod and obs on mean model depth (depth_lev)
    """
    time_lev = time_mod/86400
    depth_lev = depth_mod
    
    obs_st = gcal2jd(1760,1,1)[0]+gcal2jd(1760,1,1)[1]+0.5
    mod_st = gcal2jd(1979,9,3)[0]+gcal2jd(1979,9,3)[1]+0.5
    delta_cal = mod_st - obs_st
    time_obs_1979 = time_obs -delta_cal

    # 1.
    temp_obs_tinterp = np.zeros((len(depth_obs),len(time_lev)))
    for index in range(len(depth_obs)):
        temp_obs_tinterp[index,:] = np.interp(time_lev,\
                                          time_obs_1979,\
                                          temp_obs[index,:])

    # 2
    temp_obs_interp = np.zeros((len(depth_lev),len(time_lev)))
    temp_mod_interp = temp_mod  #[1:,:] #np.zeros((len(depth_lev),len(time_lev)))

    for index in range(len(time_lev)):
        temp_obs_interp[:,index] = np.interp(depth_lev,\
                                depth_obs,\
                                temp_obs_tinterp[:,index],\
                                left=None, right=None)

        #temp_mod_interp[:,index] = np.interp(depth_lev,\
                #                        depth_mod[:],\
                #                temp_mod[:,index])
    
    time_1950 = gcal2jd(1950,1,1)[0]+gcal2jd(1950,1,1)[1]+0.5
    delta_cal = time_1950 - mod_st
    time_lev_1950 = time_lev - delta_cal
    
    return [temp_mod_interp,temp_obs_interp,time_lev_1950,depth_lev]


def NEMO_fancy_datestr( time_counter, time_origin):
    """
    # Function: Take number input. Convert to a formatted date string referenced to a specific data.
    #    Input:
    #        time_counter - seconds since 1950
    #        time_origin - attributes for input: refence time
    #           time_origin = dataset4.variables['time_counter'].time_origin
    #    Output:
    #        time_str - fancy string for plotting 
    #        time_datetime - date as a datenumber
    #        flag_err - [1/0]: [1] if the reference date is not readable
    #           [0] if the reference date is readable, everything is OK
    #    Useage:
    #        [time_str, time_datetime, flag_err] = NEMO_fancy_datestr( time_counter, time_origin )
    """

    origin_datetime = datetime.datetime.strptime(time_origin, '%Y-%m-%d %H:%M:%S')
    time_datetime = [ origin_datetime + datetime.timedelta(seconds=time_counter[i]) for i in range(len(time_counter)) ]
    time_str =  [ datetime.datetime.strftime(time_datetime[i], '%Y-%m-%d %H:%M:%S') for i in range(len(time_counter)) ] 
    flag_err = 0 
    return [time_str, time_datetime, flag_err]

def doodsonX0(dates, elevation):
    """
#==============================================================================        
    # The Doodson X0 filter is a simple filter designed to damp out the main tidal frequencies. 
    # It takes hourly values, 19 values either side of the central one. A weighted average is taken with the following weights
    #  (1010010110201102112 0 2112011020110100101)/30.
    #  http://www.ntslf.org/files/acclaimdata/gloup/doodson_X0.html
    #
    # Learning how to apply a Doodson filter - starting from https://github.com/pwcazenave/PyFVCOM/blob/master/PyFVCOM/tappy.py
    #
    # In "Data Analaysis and Methods in Oceanography":
    #
    # "The cosine-Lanczos filter, the transform filter, and the
    # Butterworth filter are often preferred to the Godin filter,
    # to earlier Doodson filter, because of their superior ability
    # to remove tidal period variability from oceanic signals."
    """
    kern = [1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 2, 0, 1, 1, 0, 2, 1, 1, 2,
            0,
        2, 1, 1, 2, 0, 1, 1, 0, 2, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1]

    half_kern = len(kern)//2
    nslice = slice(half_kern, -half_kern)

    kern = [i/30.0 for i in kern]
    relevation = np.apply_along_axis(lambda m: np.convolve(m, kern, mode=1), axis=0, arr=elevation)

        # Pad out filter signal with NaNs
    relevation = np.lib.pad(relevation[nslice,:,:], ((19,19),(0,0),(0,0)), 'constant', constant_values=(np.nan))

    #relevation = np.convolve(elevation, kern, mode = 1)
    #return dates[nslice], relevation[nslice,:,:]
    return relevation

def getNEMObathydepth(lat0,lon0, filename='/Users/jeff/DATA/anyTide/NEMO/bathy_fix_AMM60.rmax.04.nc'):
    """
    Import a bathy depth from given coordinates
    """

    f = netcdf.netcdf_file(filename, 'r')
    # load in temperature, salinity 
    bathy = f.variables['Bathymetry'].data # (y, x) 

    #load lat and lon
    lat = f.variables['lat'].data # (y,x)
    lon = f.variables['lon'].data # (y,x)

    f.close()
    diff = abs(lat-lat0) + abs(lon-lon0)
    return bathy[np.where( diff == diff.min())]

def delta_diagnose( profile, time_counter, depth, max_depth ):
    """ 
    INPUT:
    profile data
    time_counter - seconds since 1950
    depth - metres (+ve) from the surface (z=0)
    max_depth - depth overwhich integral is calculated
    
    OUTPUT:
    pycnocline depth (m) - instantaneous
    pycnocline depth tidal filtered (m) - instantaneous
    pycnocline variance (m^2) - instantaneous
    time in python datetime speak - instantaneous
    
    pycnocline depth (m) - running window processed
    pycnocline variance (m^2) - running window processed
    time in python datetime speak - running window processed 

    Note. Initially since my running windows are actually 3 day chunks running window time will be different to instantaneous time
    When running window analysis works, these two time products should be the same.
    
    Assume the time_counter is in seconds since 1950 as is usual, not days since 1950.
    
    Usage:
    [delta, delta_nt, delta_var, time_datetime,  delta_runwin, delta_var_runwin, time_datetime_runwin] = delta_diagnose( profile, time_counter, depth, max_depth )
    """
    # Assume time counter data is relative to 1950 and in seconds
    time_origin = '1950-01-01 00:00:00'


    # Compute thermocline depth as first moment of stratification
    [nz,nt] = np.shape(profile) 
    
    # Find the index for the shallowest of 200m or bed
    index_bed = np.sum(np.isfinite(profile), axis = 0, dtype='int')-1 # Sum over depth. Returns number of finite levels per timestep
    index_200m = np.tile(np.argmin( abs(depth - max_depth)), nt)
    index = [ min([nz-2, index_200m[k], index_bed[k]]) for k in range(nt)] # nz-2 avoids dz being computed with levels that don't exist
    temp_bot = np.zeros(nt) # the temperature at min(depth,200m)
    temp_bot[:] = [profile[index[mm], mm] for mm in range(nt)]

    dz = np.zeros((nz,nt))
    temp_bar = np.zeros((nt))
    for m in range(nt):
        for k in range(index[m]):
            dz[k+1,m] = depth[k+1] - depth[k] # depth and dz are POSITIVE. k=0 is surface layer                
        temp_bar[m] = np.nansum( dz[:,m] * profile[:,m] ) / ( depth[index[m]] )

    temp_top = np.zeros(nt)
    temp_top[:] = profile[0, :]
    delta = depth[index] * (temp_bot - temp_bar) / (temp_top - temp_bot)
    delta = -np.tile(delta, (1,1,1)).T # DoodsonX0 expects a (time, y, x) array 
    
    # Expect time data in seconds
    if min( time_counter[1:] - time_counter[0:-1]) < 59: 
        print 'Time data is probably not in seconds. Expect hourly data in seconds since 1950.'
        return
    if np.nanmean( time_counter[1:] - time_counter[0:-1]) < 59: 
        print 'Time data is probably not hourly. Expect hourly data in seconds since 1950.'
        return

    #    [time_str,time_datetime, flag_err] = NEMO_fancy_datestr( np.array(time_counter[:]*24*3600, dtype='int'),time_origin ) # internal tide data
    [time_str,time_datetime, flag_err] = NEMO_fancy_datestr( np.array(time_counter[:], dtype='int'),time_origin ) # internal tide data


    
    
    # Tidally filter pycnocline depth with Doodson filter
    ################################
    # Tile the delta to increase the timeseries length if it is less than the Doodson filter length.
    if np.shape(profile)[0] < 19*2+1:
        delta = np.tile(delta, (2,1,1));
        print 'time series too short to do Doodson filter - tiled in time'

    ## Apply DoodsonX0 filter
    # For this function the time is assumed to be hourly
    delta_nt = doodsonX0(time_counter,delta);
    ## Define the internal tide variance
    ################################
    internal_tide_map = np.nanvar(delta - delta_nt, axis = 0)
    ## Define the pycnocline depth
    pycn_depth_map = np.nanmean(abs(delta_nt), axis = 0)


    ## Define the internal tide variance in 3 day chunks
    ####################################################
    [nz,nt] = np.shape(profile) # temp_mod and temp_obs are same size
    i = 0
    internal_tide_map_3day = np.zeros((int(nt/(24*3)), 1,1))
    pycn_depth_map_3day    = np.zeros((int(nt/(24*3)), 1,1))
    time_counter_3day      = np.zeros((int(nt/(24*3)), 3*24))
    time_datetime_3day  = np.array([datetime.datetime(1900,1,1) for loop in xrange(int(nt/(24*3)))]) # dummy datetime array
    #print time_datetime_3day[jj][ii]
    while ((3*i+3)*24 <= nt):
            internal_tide_map_3day[i,:,:] = np.nanvar(delta[3*i*24:(3*i+3)*24,:,:]  - delta_nt[3*i*24:(3*i+3)*24,:,:], axis = 0)
            pycn_depth_map_3day[i,:,:]    = np.nanmean( abs(delta_nt[3*i*24:(3*i+3)*24,:,:]), axis = 0)
            time_counter_3day[i,:] = time_counter[3*i*24:(3*i+3)*24]
            time_datetime_3day[i] = time_datetime[int(3*i*24 + 3*24 // 2)] # store middle times
            #print time_datetime[jj][ii][int(3*i*24 + 3*12)]
            #print time_datetime_3day[jj][ii][i]
            i += 1
#                print 'Chunking pycnocline data: ',i,' of ',int(nt/(24*3)),'. [jj,ii,nt]=',int(jj),int(ii),int(nt)

    print 'Chunking done.'
    # Do some relabelling for new notation that I want to implement in the above
    delta_var = np.squeeze(internal_tide_map)
    delta_runwin = np.squeeze(pycn_depth_map_3day)
    delta_var_runwin = np.squeeze(internal_tide_map_3day)
    time_datetime_runwin = np.squeeze(time_datetime_3day)
    delta=np.squeeze(delta)
    delta_nt=np.squeeze(delta_nt)


    return [delta, delta_nt, delta_var, time_datetime,  delta_runwin, delta_var_runwin, time_datetime_runwin]

def readMODELnc(filename, var):
    """ 
    Read a variable from a NEMO output (netcdf 3 or 4)
    """

    f = netcdf.netcdf_file(filename, 'r')
    data = f.variables[var].data
    f.close()

    return data



def appendNC(ncname,time_array,depth_array,longitude,latitude,tmod,tobs):
    """
    Appends a netcdf file with new variable
    """
    f = netcdf.netcdf_file(ncname, 'w')
    f.history = 'Temperature along FASTNEt mooring, computed from bodc files'

    f.createDimension('time',len(time_array))
    f.createDimension('depth',len(depth_array))
    f.createDimension('coord',1)

    time = f.createVariable('time', 'f', ('time',))
    time[:] = time_array
    time.units = 'days since 1950-01-01'

    depth = f.createVariable('depth','f',('depth',))
