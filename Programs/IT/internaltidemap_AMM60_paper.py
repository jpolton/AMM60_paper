# %load internaltidemap_AMM60_paper.py
#
# internaltidemap_AMM60_paper.py
#
# Diagnostics for the internal tide. Aim: make maps of amplitude and variance of IT for z<200m
#
# Origin: internaltidemap_AMM60_paper.ipynb
# jpolton 2/11/16

from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import datetime

import os # Note sure I use this
import glob # For getting file paths
import copy # For deep copying variables
from amm60_data_tools import getNEMObathydepth
from amm60_data_tools import NEMO_fancy_datestr
from amm60_data_tools import doodsonX0

import matplotlib.pyplot as plt  # plotting
%matplotlib inline


##############################################################################
# Check host name and username
import socket
hostname = socket.gethostname()

import getpass
username = getpass.getuser()

if 'livmaf' in hostname and username in ['jeff','jelt']:
    dirroot = '/Volumes'
    speedflag = True # only load in one file
else:
    dirroot = ''
    speedflag = False

##############################################################################
# Define functions for plotting multi subplot figure 
##############################################################################

# Plot the time mean depth and density fields
#############################################
def plotit(x,y,var,label,isubplot):
    plt.subplot(2, 2, isubplot)
    plt.pcolormesh(x, y, var )
    plt.colorbar()
    plt.ylabel('lat')
    plt.xlabel('lon')
    plt.title(label)


# Plot map of thermocline variance and depth
############################################
def plotit_sub(x,y,var,label,clim,s_subplot):
    plt.subplot(int(s_subplot[0]), int(s_subplot[1]), int(s_subplot[2]) )
#    plt.pcolormesh(x, y, var, cmap='nipy_spectral')
#    plt.colorbar()


    cs = plt.pcolormesh(x,y,var, cmap=plt.cm.gnuplot)
    cs.cmap.set_under('grey')
    cs.set_clim(clim)
    cb = plt.colorbar(cs, extend="max") # Extend the upper end of the colorbar    

    plt.ylabel('lat')
    plt.xlabel('lon')
    plt.ylim([+45,+63])
    plt.xlim([-14,+14])
    plt.title(label)
    




##############################################################################
# Set stuff up
##############################################################################
dirnameK = dirroot+'/projectsa/FASTNEt/kariho40/AMM60/RUNS/2010_2013/IT/' # 2010 diagnostics
dirname = dirroot+'/projectsa/FASTNEt/jelt/AMM60/RUNS/2010_2013/IT/' # 2012 diagnostics
region = 'Celtic'
#region = 'Malin'
#region = 'NSea'

data_flag = '200m'

#filename = 'AMM60_1h_20100306_20100306_diagIT_grid_T_200m.nc' # Original (March) with mask problems
#filename = 'AMM60_1h_20100306_20100307_diagIT_grid_T.nc' # New (March) with fixed mask
#filename = 'AMM60_1h_20100704_20100708_diagIT_grid_T.nc' # July 2010
#filename = 'AMM60_1h_20120601_20120610_diagIT_grid_T.nc' # June A 2012
#filename = 'AMM60_1h_20120611_20120620_diagIT_grid_T.nc' # June B 2012
#filename = 'AMM60_1h_20120621_20120630_diagIT_grid_T.nc' # June C 2012
#filename = 'AMM60_1h_20120701_20120710_diagIT_grid_T.nc' # July A 2012
#filename = 'AMM60_1h_20120711_20120720_diagIT_grid_T.nc' # July B 2012
#filename = 'AMM60_1h_20120721_20120730_diagIT_grid_T.nc' # July C 2012
#filename = 'AMM60_1h_20120731_20120809_diagIT_grid_T.nc' # July D 2012
varsurf = 'rhop_surf'
varbot = 'rhopc_bot'
varave = 'rhop_ave'
vardep = 'depth_rhopc'
varnlev = 'nlev'
print 'Diagnostics are calculated over the upper 200m'



##############################################################################
## Load in static files: bathymetry and nlev data
# the actual bathymetry (not the stuff that gets clipped at 200m)
##############################################################################
b = Dataset(dirroot+'/projectsa/FASTNEt/kariho40/AMM60/BATHY/GEBCO2014/BATHY_NOSMOOTH/bathyfile_AMM60.nc')

bathy = b.variables['Bathymetry'][:] # (y, x)
#load lat and lon
blat = b.variables['lat'][:] # (y,x)
blon = b.variables['lon'][:] # (y,x)

b.close()

## Load in the vertical levels nlev diagnosics from an old file
g = Dataset(dirnameK+'AMM60_1h_20100704_20100708_diagIT_grid_T.nc')
if data_flag == '200m':
    nlev = g.variables['nlev'][:] #(time_counter, y, x)
g.close()

print 'np.shape(blat): ',np.shape(blat)
    
    

##############################################################################
# Find all the files and loop over them
##############################################################################
#define some arrays

rho_top = []
rho_bot = []
rho_bar = []
time_counter = []
H = []
first = True

if speedflag == True:
    filenames = glob.glob(dirname+'AMM60_1h_2012????_201208??_diagIT_grid_T.nc')
    print 'fix to only read in one file on laptop'
else:
    filenames = glob.glob(dirname+'AMM60_1h_2012????_2012????_diagIT_grid_T.nc')

for fullfilepath in filenames:

    filename = fullfilepath.replace(dirname,'')
    ## Read in the data
    print 'Reading: '+filename+'\n'

    f = Dataset(fullfilepath)

    # load in pycnocline variables. NEED the [:] to slice the NetCDF data, o/w it is not manipulatable
    rhop_surf = f.variables[varsurf][:] # (time_counter, deptht, y, x)
    rhop_ave = f.variables[varave][:] # (time_counter, deptht, y, x)
    rhop_bot = f.variables[varbot][:] # (time_counter, deptht, y, x)
    #if data_flag == '200m':  # Shifted this to separate call since newer files don't have this variable
    #    nlev = f.variables[varnlev][:] #(time_counter, y, x)

    # load in time and depth
    depth = f.variables[vardep][:] # (time_counter, deptht, y, x)
    f_time_counter = f.variables['time_counter'][:] # vector
    time_origin = f.variables['time_counter'].time_origin
    time_calendar = f.variables['time_counter'].calendar
    time_units = f.variables['time_counter'].units

    #load lat and lon
    nav_lat = f.variables['nav_lat'][:] # (y,x)
    nav_lon = f.variables['nav_lon'][:] # (y,x)


    f_rho_top = np.squeeze(rhop_surf);
    f_rho_bot = np.squeeze(rhop_bot);
    f_rho_bar = np.squeeze(rhop_ave);
    f_H = np.squeeze(depth);

    f.close() # close the netcdf mooring file


    
    # Bundle all the data along the time axis
    #########################################
    if (first):
        time_counter = f_time_counter
        H            = f_H
        rho_top      = f_rho_top
        rho_bot      = f_rho_bot
        rho_bar      = f_rho_bar
        first = False
    else:
        time_counter = np.append(time_counter, f_time_counter, axis=0)
        H            = np.append(H,f_H, axis=0)
        rho_top      = np.append(rho_top, f_rho_top, axis=0)
        rho_bot      = np.append(rho_bot, f_rho_bot, axis=0)
        rho_bar      = np.append(rho_bar, f_rho_bar, axis=0)

        
        
# Load is SSH data for ST4
##########################
#fullfilepath1 = '/projectsa/FASTNEt/kariho40/AMM60/RUNS/D376/AMM60_1h_20120504_20120610_fastnet_ST4_grid_T.nc'
#fullfilepath2 = '/projectsa/FASTNEt/kariho40/AMM60/RUNS/D376/AMM60_1h_20120611_20120702_fastnet_ST4_grid_T.nc'

# Though these files are from a different simulation they overlay the common periods very well, though not perfectly.
# It is good enough for these purposes which are visual.
fullfilepath1 = dirroot+'/projectsa/FASTNEt/kariho40/AMM60/RUNS/2010_2013/NO_DIFF/AMM60_1h_20120421_20120619_fastnet_ST4_grid_T.nc'
fullfilepath2 = dirroot+'/projectsa/FASTNEt/kariho40/AMM60/RUNS/2010_2013/NO_DIFF/AMM60_1h_20120620_20120818_fastnet_ST4_grid_T.nc'


f = Dataset(fullfilepath1)

# load in time and depth
depth_ST4 = f.variables['depth'][:] # (time_counter, deptht, y, x)
time_counter_ST4  = f.variables['time_counter'][:]

f = Dataset(fullfilepath2)

# load in time and depth
depth_ST4 = np.append(depth_ST4,f.variables['depth'][:] , axis = 0) # (time_counter, deptht, y, x)
time_counter_ST4  = np.append( time_counter_ST4, f.variables['time_counter'][:], axis = 0)

# Remove the mean from the depth
depth_ST4 = depth_ST4 - np.mean(depth_ST4, axis = 0)


##############################################################################
# Process data
##############################################################################
print 'Process data'
print '############'

# Process the time data
################################
# Note that the Error flag doesn't work and I haven't actually checked it. What happens with leap years etc...
[time_str,     time_datetime,     flag_err] = NEMO_fancy_datestr( time_counter,     time_origin ) # internal tide data
[time_str_ST4, time_datetime_ST4, flag_err] = NEMO_fancy_datestr( time_counter_ST4, time_origin ) # SSH data



## Compute delta
################################
delta={};
delta = H*( rho_bot - rho_bar) / (rho_top - rho_bot)

## Compute stratification
################################
strat = {};
strat = (rho_top - rho_bot) / H;


# Define masks
###############################
import numpy.ma as ma

[nt,ny,nx] = np.shape(strat)
mask_shelf = np.reshape( ma.masked_where(nlev[0,:,:] < 41, nlev[0,:,:]) ,(ny,nx) )

mask_land = np.reshape( ma.masked_where(nlev[0,:,:] == 0, np.ones((ny,nx))) ,(ny,nx) )

mask_200m = np.reshape( ma.masked_where(bathy[:,:] >= 200, np.ones((ny,nx))) ,(ny,nx) )

# Mask NW corner and west of 12W
m = (63.-56.)/(0.-(-14.)) # Gradient of corner slice
c = 63. # Slice lat intercept when lon=0
t2 = np.array(nav_lat - m*nav_lon - c >= 0 ,dtype=bool) # Mask NW corner slice
t1 = np.array(nav_lon < -12. , dtype=bool) # Mask west of 12W
mask_corners = np.reshape( ma.masked_where( t1+t2 , np.ones((ny,nx))) ,(ny,nx) ) # Boolean addition is OR not SUM

# Agregate masks associated with geographic domain of interest. Exclude strat mask.
mask = mask_land*mask_200m*mask_corners # Ones and zeros, wet domain mask


mean_strat = np.mean( strat ,axis=0)
mask_strat = ( mean_strat >= -2E-3 ).astype(int)*(-999)   # Good vals: 0 / bad vals: -999. Plus mask_strat to variable when plotting.
#mask_strat = ma.masked_where(np.mean( strat ,axis=0) >= -2E-3, np.ones((ny,nx)))



# Tidally filter pycnocline depth with Doodson filter
################################
# Tile the delta to increase the timeseries length if it is less than the Doodson filter length.
if np.shape(delta)[0] < 19*2+1:
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
## Define a potential energy of the pycnocline disturbance
#pe = 0.5 * 9.81 * np.mean((rho_top - rho_bed) * ( internal_tide_map**2 / H ), axis=0)

#for (x, y), element in np.ndenumerate(nav_lat):
#    print(x, y, internal_tide_map[x,y], pycn_depth_map[x,y], np.mean(H[:,x,y],axis=0), nav_lat[x,y],nav_lon[x,y])


## Define the internal tide variance in 3 day chunks
####################################################
i = 0
internal_tide_map_3day = np.zeros((nt/(24*3), ny,nx))
pycn_depth_map_3day    = np.zeros((nt/(24*3), ny,nx))
time_counter_3day      = np.zeros((nt/(24*3), 3*24))
while ((3*i+3)*24 < nt): 
        internal_tide_map_3day[i,:,:] = np.nanvar(delta[3*i*24:(3*i+3)*24,:,:]  - delta_nt[3*i*24:(3*i+3)*24,:,:], axis = 0)
        pycn_depth_map_3day[i,:,:]    = np.nanmean( abs(delta_nt[3*i*24:(3*i+3)*24,:,:]), axis = 0)
        time_counter_3day[i,:] = time_counter[3*i*24:(3*i+3)*24]
        i += 1
        print 'Chunking pycnocline data: ',i,' of ',(nt/24)/3







##############################################################################
# Plot density and depth fields
##############################################################################
print 'Plot stuff'
print '##########'

fig, ax = plt.subplots(2,2)
plt.rcParams['figure.figsize'] = (20.0, 16.0)


## analysis depth
###############################
var = np.log10(np.mean(H,axis=0))
plotit(nav_lon,nav_lat,var,'log10(water column analaysis depth)',1)

## top density
###############################
var = np.mean(rho_top,axis=0)
var[var==0]=np.nan
clim = [np.nanpercentile(var, 5), np.nanpercentile(var, 95)]
#print 'colorbar limits:', clim
#var = ma.masked_where(var <= 1024, var)
plotit(nav_lon,nav_lat,var,'surface density',2)
plt.clim(clim)

## bottom density
###############################
var = np.mean(rho_bot,axis=0)
var[var==0]=np.nan
#var = ma.masked_where(var <= 1020, var)
plotit(nav_lon,nav_lat,var,'density at deepest analysis depth',3)
plt.clim(clim)

## mean density, rho_bar
###############################
#var2=( rho_bot - rho_bar) / (rho_top - rho_bot)*np.mean(H,axis=0)
#var = var2[0,:,:]
var=np.mean(rho_bar,axis=0)
var[var==0]=np.nan
clim = [np.nanpercentile(var, 5), np.nanpercentile(var, 95)]
#var = ma.masked_where(var <= 1024, var)
plotit(nav_lon,nav_lat,var,'ave density over analysis depth',4)
plt.clim(clim)


##############################################################################
# Plot pycnocline statistics
##############################################################################

fig, ax = plt.subplots(2,2)
plt.rcParams['figure.figsize'] = (20.0, 16.0)
plt.suptitle(filename)

## Analysis depth range
###############################
var = np.mean(H,axis=0)*mask_land*mask_200m
clim = [np.nanpercentile(var, 5), np.nanpercentile(var, 95)]
#print 'subplot 1: percentile range:',clim
plotit_sub(nav_lon,nav_lat,var,'Analysis depth range (m)',clim,'221')
#plt.clim(clim)
#plt.contour(nav_lon,nav_lat,depth[0,:,:]*mask_shelf, [100,200])

## pycnocline depth
###############################
#var = pycn_depth_map
var = delta[0,:,:]*mask_land*mask_200m # *mask_strat
var[var>=200]=np.nan
clim = [-50, -10]
#clim = [np.nanpercentile(var, 5), np.nanpercentile(var, 95)]
#print 'subplot 2: percentile range:',clim
plotit_sub(nav_lon,nav_lat,var+mask_strat,'pycnocline depth (m)',clim,'222')
#plt.clim(clim)
#plt.contour(nav_lon,nav_lat,H[0,:,:]*mask_shelf, [100,200])


## average bulk stratification
###############################
var = np.mean( -strat ,axis=0)*mask_land*mask_200m #*(mask_strat*mask_strat)
#var = np.max(nlev,axis=0) - np.min(nlev,axis=0)
clim = [0, 0.02]
#clim = [0, np.nanpercentile(var, 95)]
#print 'subplot 4: percentile range:',clim
plotit_sub(nav_lon,nav_lat,var+mask_strat,'-mean bulk stratification (kg/m^4)',clim,'223')
#plt.clim(clim)


## pycnocline depth variance
###############################
var = np.log10(internal_tide_map*mask_land*mask_200m) #*mask_strat)
clim = [0., 2.6]
#clim = [np.nanpercentile(var, 1), np.nanpercentile(var, 99)]
#print 'subplot 3: percentile range:',clim
plotit_sub(nav_lon,nav_lat,var+mask_strat,'log10[pycnocline depth variance (m)]',clim,'224')
#plt.clim(clim)





##############################################################################
# Plot pycnocline statistics in 3 day chunks
##############################################################################
#for i in range(6,7):
for i in range(nt/(24*3)):

    # Find indices in SSH ST4 data that correspond to the IT data
    ind = [ii for ii in range(len(time_counter_ST4)) if time_counter_ST4[ii] in time_counter_3day[i,:]]

    fig, ax = plt.subplots(2,2)
    plt.rcParams['figure.figsize'] = (20.0, 16.0)
    plt.suptitle(str(time_datetime_ST4[ind[0]])+' to '+str(time_datetime_ST4[ind[-1]]), fontsize=20)

      ## SSH
    ###############################
    var = depth_ST4
    hlim = [-0.25, 0.25]
    plt.subplot(4,1,1)
    plt.plot(time_datetime_ST4, depth_ST4[:,-1,1,1])
    plt.plot([time_datetime_ST4[ii] for ii in ind], depth_ST4[ind,-1,1,1], 'r')
    dstart = datetime.datetime(2012,6,1)
    dend = datetime.datetime(2012,8,9)
    plt.xlim(dstart, dend)
    plt.xlabel('time')
    plt.ylabel('SSH above mean (m)')
    plt.title('ST4 SSH')
    #plt.clim(hlim)

      
    ## pycnocline depth
    ###############################
    var = copy.deepcopy(pycn_depth_map_3day[i,:,:])
    clim = [0, 50]
    #clim = copy.deepcopy([0, np.nanpercentile(var*mask, 95)])
    var[mask*mask_strat==-999] = -999 # These will be off the bottom of the colorbar scale, assign grey 
    plotit_sub(nav_lon,nav_lat,var*mask,'pycnocline depth (m)',clim,'323')

    ## pycnocline depth histogram
    ###############################
    plt.subplot(3,2,5)
    plt.plot( np.sort(pycn_depth_map_3day[i,:,:]*mask, axis=None) )
    plt.xlabel('sorted grid box index')
    plt.ylabel('pycnocline depth (m)')
    plt.ylim(clim)
    plt.xlim(0,3E5)

    
    
    ## pycnocline depth variance
    ###############################    
    var = copy.deepcopy(internal_tide_map_3day[i,:,:])
    #clim = [0, 30]
    clim = [np.log10(1), np.log10(30)]
    #clim = copy.deepcopy([0, np.nanpercentile(var*mask, 95)])
    var[mask*mask_strat==-999] = -999 # These will be off the bottom of the colorbar scale, assign grey 
    #plotit_sub(nav_lon,nav_lat,var*mask,'pycnocline depth variance(m)',clim,'322')
    plotit_sub(nav_lon,nav_lat,np.log10(var)*mask,'log10[pycnocline depth variance(m)]',clim,'324')
        
    ## pycnocline depth variance histogram
    ###############################
    plt.subplot(3,2,6)
    #plt.plot( np.sort(internal_tide_map_3day[i,:,:]*mask, axis=None) )
    #plt.ylabel('pycnocline depth variance (m)')  
    plt.plot( np.sort(np.log10(internal_tide_map_3day[i,:,:])*mask, axis=None) )
    plt.ylabel('log10[pycnocline depth variance (m)]')  
    plt.xlabel('sorted grid box index')
    plt.ylim(clim)
    plt.xlim(0,3E5)
        


    ## Save output
    ###############################
    fname = dirroot+'/scratch/jelt/tmp/internaltidemap'+str(i).zfill(3)+'_IT_stats.png'
    plt.savefig(fname)