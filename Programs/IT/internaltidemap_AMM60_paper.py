# %load internaltidemap_AMM60_paper.py
#
# internaltidemap_AMM60_paper.py
#
# Diagnostics for the internal tide. Aim: make maps of amplitude and variance of IT for z<200m
#
# Origin: internaltidemap_AMM60_paper.ipynb
# jpolton 2/11/16
if(1): # Change flag to skip data loading and processing step if already done.
    from netCDF4 import Dataset
    import numpy as np
    import numpy.ma as ma # masks
    import datetime

    import os # Note sure I use this
    import glob # For getting file paths
    import copy # For deep copying variables
    
    import sys
    sys.path.append('../FASTNEt/') # Add the directory with the amm60_data_tools.py file to path

    from amm60_data_tools import NEMO_fancy_datestr # convert NEMO time 
    from amm60_data_tools import delta_diagnose # compute pycnocline depth and variance. Incl running window filter
    from amm60_data_tools import window_strat # compute running window filtered stratification

    import matplotlib.colors as colors # colorbar with log distribution
    import matplotlib.pyplot as plt  # plotting
    %matplotlib inline


    ##############################################################################
    # Check host name and username.
    # Depending on machine and modify path tree and flag to load partial data set
    import socket
    hostname = socket.gethostname()

    import getpass
    username = getpass.getuser()

    if 'livmaf' in hostname and username in ['jeff','jelt']:
        dirroot = '/Volumes'
        speedflag = True # only load in one file - speeds things up when debugging
    elif 'livljobs' in hostname and username in ['jeff','jelt']:
        dirroot = ''
        speedflag = False # only load in one file
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

        cs = plt.pcolormesh(x,y,var, cmap=plt.cm.gnuplot)
        cs.cmap.set_under('grey')
        cs.set_clim(clim)
        cb = plt.colorbar(cs, extend="max") # Extend the upper end of the colorbar    

        plt.ylabel('lat')
        plt.xlabel('lon')
        plt.ylim([+45,+63])
        plt.xlim([-14,+14])
        plt.title(label)


    # Plot map of thermocline variance and depth. (As above but on log colour scale)
    ############################################
    def plotit_sub_log(x,y,var,label,clim,s_subplot):
        plt.subplot(int(s_subplot[0]), int(s_subplot[1]), int(s_subplot[2]) )

        cs = plt.pcolormesh(x,y,var, cmap=plt.cm.gnuplot,
                            norm=colors.LogNorm(vmin=clim[0], vmax=clim[1]))
        cs.cmap.set_under('grey')
        cb = plt.colorbar(cs, extend="max") # Extend the upper end of the colorbar    
        cb.set_ticks(range(1,11))
        cb.set_ticklabels(range(1,11))

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

    varsurf = 'rhop_surf'
    varbot = 'rhopc_bot'
    varave = 'rhop_ave'
    vardep = 'depth_rhopc'
    varnlev = 'nlev'

    data_flag = '200m' # loaded deep data is outputted at 200m or bed.
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

    del f_rho_top
    del f_rho_bar
    del f_rho_bot
    del f_time_counter
    del f_H

    # Load is SSH data for ST4
    ##########################
    # NB Though these files are from a different simulation they overlay the common periods very well
    # It is good enough for these purposes which are visual, to indicate the spring neap phase
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

        
    # Stack rho data into a 3 layer in z-direction, so data can be treated as profiles
    profile = rho_top[:,np.newaxis,:,:]
    profile = np.append(profile, rho_bar[:,np.newaxis,:,:], axis=1)
    profile = np.append(profile, rho_bot[:,np.newaxis,:,:], axis=1)

    # compute time-means then clear large variables
    mean_rho_top = np.mean(rho_top, axis=0)
    mean_rho_bar = np.mean(rho_bar, axis=0)
    mean_rho_bot = np.mean(rho_bot, axis=0)
    del rho_top
    del rho_bar
    del rho_bot
    
    ## Compute delta and variance properties with 3 day windowing
    #################################
    print 'compute pycnocline fields - This is super slow. Make some tea.'
    max_depth = 0 # not used for 3D profiles
    [delta, delta_nt, internal_tide_map, time_datetime,  pycn_depth_map_3day, internal_tide_map_3day,
                 time_datetime_3day, time_counter_3day] = delta_diagnose( profile, time_counter, H, max_depth )


    ## Compute stratification with 3 day windowing
    #################################
    print 'Time window stratification - quite slow'
    strat_3day = window_strat( profile, time_counter, H )


    profile = [] # Clear some memory

    # compute time-mean for plotting then clear space
    mean_H = np.mean(H, axis=0)
    del H

    # Process the SSH time data
    ################################
    # Note that the Error flag doesn't work and I haven't actually checked it. What happens with leap years etc...
    [time_str_ST4, time_datetime_ST4, flag_err] = NEMO_fancy_datestr( time_counter_ST4, time_origin ) # SSH data


    # Define masks
    ###############################
    [nt,ny,nx] = np.shape(delta)
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

    [runwin_nt,ny,nx] =  np.shape(internal_tide_map_3day[:,:,:]) # ny=nx=1 for 1d profiles

    mask_strat = ( strat_3day >= -3E-3 ).astype(int)*(-9999)   # Good vals: 0 / bad vals: -9999. 
    mean_strat = np.mean( strat_3day ,axis=0)
    mean_mask_strat = ( mean_strat >= -3E-3 ).astype(int)*(-9999)   # Good vals: 0 / bad vals: -9999.

    # Process sorted variances
    ###############################
    print 'process sorted variances'

    # Define new array to store sorted variance data
    sortvar = np.zeros((runwin_nt,100)) # array [nt,% of finite area] of domain with var at value. 
    for i in range(runwin_nt):
        # Sort snapshot of data 
        var = mask*copy.deepcopy(internal_tide_map_3day[i,:,:]) # function of [t=0,y,x] * (spatial mask)
        var[mask*mask_strat[i,:,:]==-9999] = np.nan # nan for unstratified locations
        tt = np.sort(np.log10(var), axis=None)

        # Interpolate onto 100 points, excluding nans
        not_nan = ~np.isnan(tt)
        indices = np.arange(len(tt[not_nan]))
        ind_short = np.linspace(0,np.sum(not_nan)-1,100)
        sortvar[i,:] = np.interp(ind_short, indices, tt[not_nan])


#""""""
#End of comment to skip data loading and processing
#""""""

##############################################################################
# Plot density and depth fields
##############################################################################
print 'Plot stuff'
print '##########'

fig, ax = plt.subplots(2,2)
plt.rcParams['figure.figsize'] = (20.0, 16.0)


## analysis depth
###############################
var = np.log10(mean_H)
plotit(nav_lon,nav_lat,var,'log10(water column analaysis depth)',1)

## top density
###############################
var = mean_rho_top
var[var==0]=np.nan
clim = [np.nanpercentile(var, 5), np.nanpercentile(var, 95)]
plotit(nav_lon,nav_lat,var,'surface density',2)
plt.clim(clim)

## bottom density
###############################
var = mean_rho_bot
var[var==0]=np.nan
plotit(nav_lon,nav_lat,var,'density at deepest analysis depth',3)
plt.clim(clim)

## mean density, rho_bar
###############################
var=mean_rho_bar
var[var==0]=np.nan
clim = [np.nanpercentile(var, 5), np.nanpercentile(var, 95)]
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
var = mean_H*mask
clim = [np.nanpercentile(var, 5), np.nanpercentile(var, 95)]
plotit_sub(nav_lon,nav_lat,var,'Analysis depth range (m)',clim,'221')

## pycnocline depth
###############################
var = delta[0,:,:]*mask
var[var>=200]=np.nan
clim = [10, 50]
plotit_sub(nav_lon,nav_lat,var+mean_mask_strat,'pycnocline depth (m)',clim,'222')


## average bulk stratification
###############################
var = np.mean( np.abs(strat_3day) ,axis=0)*mask
clim = [3E-3, 0.02]
plotit_sub(nav_lon,nav_lat,var+mean_mask_strat,'-mean bulk stratification (kg/m^4)',clim,'223')

## pycnocline depth tidal std
###############################
#var = np.log10(internal_tide_map*mask_land*mask_200m) 
var = np.sqrt(internal_tide_map)*mask
clim = [0., 10.]

if(1): # Perahps plot std with a log colourscale     
    ax3 = fig.add_subplot(2,2,4)
    cs = ax3.pcolormesh(nav_lon,nav_lat,var+mean_mask_strat, cmap=plt.cm.gnuplot,
                        norm=colors.LogNorm(vmin=0.1, vmax=10))
    cs.cmap.set_under('grey')
    cb = plt.colorbar(cs, extend="both") # Extend the upper end of the colorbar    
    ax3.set_title('pycnocline depth tidal std (m)]')
    ax3.set_ylabel('lat')
    ax3.set_xlabel('lon')
    ax3.set_ylim([+45,+63])
    ax3.set_xlim([-14,+14])

##############################################################################
# Plot pycnocline statistics in 3 day chunks
##############################################################################
#for i in range(6,7): # Demo output
for i in range(np.size(time_counter_3day[:,:][0])):

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
      
    ## pycnocline depth
    ###############################
    var = copy.deepcopy(pycn_depth_map_3day[i,:,:])
    clim = [0, 50]
    var[mask*mask_strat[i,:,:]==-9999] = -9999 # These will be off the bottom of the colorbar scale, assign grey 
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
    clim = [0.1, 10]
    var[mask*mask_strat[i,:,:]==-9999] = -9999 # These will be off the bottom of the colorbar scale, assign grey 

    var = np.sqrt(internal_tide_map_3day[i,:,:])*mask
    clim = [0., 10.]
    plotit_sub(nav_lon,nav_lat,var+mask_strat[i,:,:],'pycnocline depth tidal std (m)]',clim,'324')


    ## pycnocline depth variance histogram
    ###############################
    plt.subplot(3,2,6) 
    plt.plot( np.sort(0.5*np.log10(internal_tide_map_3day[i,:,:])*mask, axis=None) )
    plt.ylabel('log10[pycnocline depth tidal std (m)]')  
    plt.xlabel('sorted grid box index')
    plt.ylim([np.log10(0.1),np.log10(10)])
    plt.xlim(0,3E5)
        


    ## Save output
    ###############################
    fname = dirroot+'/scratch/jelt/tmp/internaltidemap'+str(i).zfill(3)+'_IT_stats.png'
    plt.savefig(fname)
    
    
    
##############################################################################    
# Plot std histogram vs time
##############################################################################

"""
Actually want to plot std(delta) rather than the variance.
Since the quantity carried is sortvar=log10(variance) need to do a little bit of algebra
variance = std^2 = 10^(sortvar)
Therefore
std = sqrt(10^sortvar) OR 10^(sortvar/2)
"""

std = np.power(10,sortvar/2)

fig = plt.figure(figsize=(12,12))

## SSH
###############################
var = depth_ST4
hlim = [-0.25, 0.25]
ax1 = fig.add_subplot(411)
ax1.plot(time_datetime_ST4, depth_ST4[:,-1,1,1])
dstart = datetime.datetime(2012,6,1)
dend = datetime.datetime(2012,8,9)
ax1.set_xlim(dstart, dend)
ax1.set_ylabel('SSH above mean (m)')
# text label
start = ax1.get_xlim()[0] + 0.5
ax1.text(start, -0.025, 'a) SSH at ST4')


## std vs time
###############################
ax2 = fig.add_subplot(412)
msh = ax2.contourf(time_datetime_3day,np.arange(100), std.T, [ 0.1, 0.5, 1,2,4,8,10],
                     cmap=plt.cm.gnuplot,
                     norm=colors.LogNorm(vmin=.1, vmax=10))
ax2.set_ylabel('domain coverage %')
ax2.set_xlabel('time')
ax2.set_xlim(dstart,dend)

# Now adding the colorbar
cbaxes = fig.add_axes([0.91, 0.125, 0.03, 0.775]) # [left, bottom, width, height]
cb = fig.colorbar(msh, cax = cbaxes, orientation='vertical', extend="both") 
# Fiddle with the colorbar ticks
cb.set_ticks([ 0.5, 1,2,4,8,10])
cb.set_ticklabels([ 0.5, 1,2,4,8,10])

# text label
start = ax2.get_xlim()[0] + 2.5
ax2.text(start, 7, 'b) std($\delta$) (m)',color='w')



# Add spring and neap snapshot std maps
#######################################

for i in [0,1]:
    count = [260, 324]
    label = ['c) ','d) ']
    # Find indices in SSH ST4 data that correspond to the IT data
    ind = [ii for ii in range(len(time_counter_ST4)) if time_counter_ST4[ii] in time_counter_3day[count[i],:]]
    ax1.plot([time_datetime_ST4[ii] for ii in ind], depth_ST4[ind,-1,1,1], 'r')


    ## pycnocline depth tidal std
    ############################### 
    std = np.sqrt(copy.deepcopy(internal_tide_map_3day[count[i],:,:]))*mask
    std[mask*mask_strat[i,:,:]==-9999] = 0.01 # These will be off the bottom of the colorbar scale, assign grey 
    
    ax3 = fig.add_subplot(2,2,3+i)
    cs = ax3.pcolormesh(nav_lon,nav_lat,std, cmap=plt.cm.gnuplot,
                        norm=colors.LogNorm(vmin=.1, vmax=10))
    cs.cmap.set_under('grey')

    ax3.set_ylabel('lat')
    ax3.set_xlabel('lon')
    ax3.set_ylim([+45,+63])
    ax3.set_xlim([-14,+14])

    # Process datetime for analysis range
    [time_str, time_datetime, flag_err] = NEMO_fancy_datestr( time_counter_3day[count[i],:], time_origin )
    ax3.text(-13, 62, label[i]
             +str(time_datetime[0].strftime('%d')+'-'+time_datetime[-1].strftime('%d %b %Y'))
             +', std($\delta$) (m)')

    
# Save output
###############################
#fname = dirroot+'/scratch/jelt/tmp/internaltidemap_std.png'
fname = 'internaltidemap_std.png'
plt.savefig(fname)