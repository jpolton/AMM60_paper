# %load fft_paper.py
#
# fft_paper.py
#
# Compute FFTs of pycnocline depth for manuscript
#
# Origin: fft_paper.ipynb
# kguihou, jpolton 15/12/16

from scipy.io import netcdf
import scipy
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
from __future__ import division
#import AMM60_tools  # I think that we have done away with this file and called it amm60_data_tools.py
%matplotlib inline

import sys
sys.path.append('../FASTNEt/') # Add the directory with the amm60_data_tools.py file to path
import amm60_data_tools # The new AMM60_tools

#########
## Set up
#########

len_seg = 192  # FFT data chunk length (data is hourly)

# Check username.
# Depending on username modify path tree
########################################

import getpass
username = getpass.getuser()

if username in ['jeff','jelt']:
    dirpath = '../FASTNEt/'
elif username in ['karen','kariho40']:
    dirpath = '/media/data1/AMM60/FASTNEt/'
else:
    dirpath = '../FASTNEt/'

print 'Karen, Hopefully your username is  correct here and your paths will be maintained. Though hopefully this should also work if you download the whole project'


## Define functions to extract the thermocline depth, compute FFT
#################################################################

def find_segment(filename):
    # Get the min/max index of available data.
    delta_obs = amm60_data_tools.readMODELnc(filename,'delta_nt_obs')
    delta_obs = np.ma.masked_where(delta_obs == 0, delta_obs)
    delta_obs = np.ma.masked_invalid(delta_obs)
    Index_start = np.argwhere(delta_obs)[0]
    Index_end = np.argwhere(delta_obs)[-1]
    
    return[Index_start,Index_end]

def compute_fft(filename,Index_start,Index_end):
    ### Load the data ###
    delta_obs = amm60_data_tools.readMODELnc(filename,'delta_obs')[Index_start:Index_end]
    delta_mod = amm60_data_tools.readMODELnc(filename,'delta_mod')[Index_start:Index_end]
    delta_3d_obs = amm60_data_tools.readMODELnc(filename,'delta_nt_obs')[Index_start:Index_end]
    delta_3d_mod = amm60_data_tools.readMODELnc(filename,'delta_nt_mod')[Index_start:Index_end]

    ### FFT ###
    print(len(delta_mod))
    fftobs = np.zeros((int(len_seg/2)+1,len(delta_mod)))
    fftmod = np.zeros((int(len_seg/2)+1,len(delta_mod)))
# Old FFT calls
#    f, fftobs = scipy.signal.welch(delta_obs-delta_3d_obs, fs=1.0, window='hanning', nperseg=len_seg, noverlap=len_seg*2/3, nfft=2*len_seg, detrend='linear', return_onesided=True, scaling='density')    
#    f, fftmod = scipy.signal.welch(delta_mod-delta_3d_mod, fs=1.0, window='hanning', nperseg=len_seg, noverlap=len_seg*2/3, nfft=2*len_seg, detrend='linear', return_onesided=True, scaling='density')
    f, fftobs = scipy.signal.welch(delta_obs-delta_3d_obs, fs=1.0, window='hanning', nperseg=len_seg, noverlap = len_seg*0.5, nfft=1*len_seg,detrend='linear', return_onesided=True, scaling='spectrum')    
    f, fftmod = scipy.signal.welch(delta_mod-delta_3d_mod, fs=1.0, window='hanning', nperseg=len_seg, noverlap = len_seg*0.5, nfft=1*len_seg,detrend='linear', return_onesided=True, scaling='spectrum')

    return[f,fftobs,fftmod]

def compute_interial_freq(filename):
    lat = amm60_data_tools.readMODELnc(filename,'lat')
    lon = amm60_data_tools.readMODELnc(filename,'lon')
    if lat < 40 and lon>40:
        print 'I think that the lat and lon are the wrong way around'
        print 'lat=',lat
        print 'lon=',lon
        lat = lon
        
        return 24/(2*np.sin( np.pi/180.*lat)) # Inertial period (hrs)
    
######################################################
## Compute the FFTs for each mooring and configuration
######################################################
    
# Loop over moorings
for mooring_num in ['ST1','ST2','ST4','ST5']:


    ###### AMM60 #####
    print('AMM60')
    filename = dirpath+'AMM60/mooring_'+mooring_num+'.nc'
    [Index_start,Index_end] = find_segment(filename)
    [famm60,fftobsamm60,fftmodamm60 ] = compute_fft(filename,Index_start,Index_end)

    print(fftobsamm60.mean())
    print(fftmodamm60.mean())



    ###### AMM7 #####
    # read first time to extract the same period as AMM60
    filename = dirpath+'AMM7/mooring_'+mooring_num+'.nc'
    [Index_start,Index_end] = find_segment(filename)
    [famm7,fftobsamm7,fftmodamm7 ] = compute_fft(filename,Index_start,Index_end)

    print(fftobsamm7.mean())
    print(fftmodamm7.mean())



    ###### NNA #####
    # read first time to extract the same period as AMM60
    filename = dirpath+'NNA/mooring_'+mooring_num+'.nc'
    [Index_start,Index_end] = find_segment(filename)
    [fnna,fftobsnna,fftmodnna ] = compute_fft(filename,Index_start,Index_end)

    print(fftobsnna.mean())
    print(fftmodnna.mean())


    ######################################################
    ## Plot data
    ######################################################

    # Plot linear y-axis
    axes = plt.figure(figsize=(25, 10),dpi=250,facecolor='w', edgecolor='k').add_subplot(111)
    #plt.figure(figsize=(25, 10), dpi=250, facecolor='w', edgecolor='k')
    plt.rcParams.update({'font.size': 25})

    period_hrs = 12.42; lineformat = 'r--'; label='M2'
    plt.plot( [1/period_hrs, 1/period_hrs], [1E-2, 1E7], lineformat)
    #ax.text(1/(3600*period_hrs), 1E6, label)

    plt.plot(famm60,fftobsamm60,color='k',label='Observations',linewidth=3)
    plt.plot(famm60,fftmodamm60,color='r',label='AMM60',linewidth=2)
    plt.plot(famm7,fftmodamm7,color='b',label='AMM7',linewidth=2)
    plt.plot(fnna,fftmodnna,color='g',label='NNA',linewidth=2)

    #period_hrs = 312.00; lineformat = 'k--'; label='M2-S2'
    #plt.semilogy( [1/period_hrs, 1/period_hrs], [1E-6, 1E7], lineformat,label='M2-S2')

    period_hrs = 26.87; lineformat = 'y--'; label='Q1'
    plt.plot( [1/period_hrs, 1/period_hrs], [1E-6, 1E7], lineformat,label='Q1')

    period_hrs = 25.82; lineformat = 'b--'; label='O1'
    plt.plot( [1/period_hrs, 1/period_hrs], [1E-6, 1E7], lineformat,label='O1')

    period_hrs = 23.93; lineformat = 'g--'; label='K1'
    plt.plot( [1/period_hrs, 1/period_hrs], [1E-6, 1E7], lineformat,label='K1')

    period_hrs = 24.07; lineformat = 'm--'; label='P1'
    #period_hrs = 14.96; lineformat = 'y--'; label='P1' # This is the phase speed in deg/hr
    plt.plot( [1/period_hrs, 1/period_hrs], [1E-6, 1E7], lineformat,label='P1')

    period_hrs = 12.66; lineformat = 'c--'; label='N2'
    plt.plot( [1/period_hrs, 1/period_hrs], [1E-6, 1E7], lineformat,label='N2')

    period_hrs = 12.42; lineformat = 'r--'; label='M2'
    plt.plot( [1/period_hrs, 1/period_hrs], [1E-6, 1E7], lineformat,label='M2')

    period_hrs = 12.00; lineformat = 'm--'; label='S2'
    plt.plot( [1/period_hrs, 1/period_hrs], [1E-6, 1E7], lineformat,label='S2')

    period_hrs = 11.967; lineformat = 'b--'; label='K2'
    plt.plot( [1/period_hrs, 1/period_hrs], [1E-6, 1E7], lineformat,label='K2')

    period_hrs = 6.21; lineformat = 'k--'; label='M4'
    plt.plot( [1/period_hrs, 1/period_hrs], [1E-6, 1E7], lineformat,label='M4')

    period_hrs = compute_interial_freq(filename); lineformat = 'k-.'; label='In' 
    plt.plot( [1/period_hrs, 1/period_hrs], [1E-6, 1E7], lineformat,label='In')


    plt.xlim([0,0.2])
    plt.ylim([0, 30])

    # Switch the freq axis ticks and labels into hours
    axes.set_xticks(1./np.arange(24,5,-1))

    axes.minorticks_off()
    axes.tick_params('x', direction='out',length=5, width=2, which='major', top='off')
    axes.tick_params('x', length=10, width=1, which='minor')

    axes.set_xticklabels(['24','','','','','','','','','','','','12','','','','','','6'])

    plt.title(mooring_num+' mooring')
    if (mooring_num == 'ST1'):
        plt.legend(fontsize=20, loc='upper right', ncol=4)
    plt.xlabel('period [hours]')
    #plt.xlabel('frequency [cyles per hour]')
    plt.ylabel('Power spectra [$m^2$]')

    name = 'fft_'+mooring_num+'.png'
    plt.savefig(str(name))

    plt.show()

