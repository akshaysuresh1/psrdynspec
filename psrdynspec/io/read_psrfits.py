# Load data from list of PSRFITS files.

from psrdynspec.io.read_header import Header
from astropy.io import fits
import numpy as np
import sys, glob
#################################################
# Extract all data from a list of PSRFITS files.
'''
Inputs:
file_list = List of .fits files
pol = List of polarization indices to extract (default = [0])
'''
def load_psrfits_data(file_list,pol=[0]):
    if len(file_list)==0:
        print('PSRFITS file list is empty. Quitting.')
        sys.exit(1)
    for i in range(len(file_list)):
        file_name = file_list[i].split('/')[-1]
        print('Reading in data from %s'% (file_name))
        data_file = fits.open(file_list[i], ignore_missing_end=True)[-1].data['DATA'][:,:,:,:,0]
        # Data shape = (Nrows, NSBLK, NPOL, NCHANS)
        data_file = data_file.reshape((data_file.shape[0]*data_file.shape[1],data_file.shape[2],data_file.shape[3]))
        # Data shape = (Ntsamples, NPOL, NCHANS)
        data_file = np.moveaxis(data_file,0,-1) # Move time samples to last axis of array.
        # Select required polarizations.
        data_file = data_file[pol]
        if (i==0):
            data = np.copy(data_file)
        else:
            data = np.concatenate((data,data_file),axis=-1)
    print('Data load complete.')
    return data
#################################################
# Read header and extract a chunk of data from a set of PSRFITS files.
'''
glob_psrfits = Glob string to parse .fits files
start_time = Start time (s) of data to load
stop_time = Stop time (s) of data to load
pol = An integer or a list of polarization indices to extract (default = [0])
'''
def extract_psrfits_datachunk(glob_psrfits, start_time, stop_time, pol=[0]):
    file_list = sorted(glob.glob(glob_psrfits))
    if len(file_list)==0:
        print('PSRFITS file list is empty. Quitting.')
        sys.exit(1)
    hdr = Header(glob_psrfits, file_type='psrfits')
    tot_time_samples = hdr.ntsamples # Total no. of time samples summed over all data in input list.
    t_samp = hdr.t_samp # Sampling time (s)

    ntsamples_per_file = np.ones(len(file_list))*hdr.subint['NSBLK']*hdr.subint['NAXIS2'] # No. of time samples per file
    # Update last entry to allow for possibility for fewer time samples in the final data set of the input list.
    # If N = len(file_list), tot_time_samples = (N-1)*ntsamples_per_file[0] + ntsamples_per_file[-1]
    ntsamples_per_file[-1] = tot_time_samples - hdr.subint['NSBLK']*hdr.subint['NAXIS2']*(len(file_list)-1)
    ntsamples_per_file = ntsamples_per_file.astype(int)
    # Cumulative no. of time samples scanned at the start of each PSRFITS file in the list.
    cumsamples_per_file = np.cumsum(ntsamples_per_file).astype(int)

    print('Supplied time range:')
    print('Start time = %.2f s'% (start_time))
    print('Stop time = %.2f s'% (stop_time))
    t_start = np.floor(start_time/t_samp).astype(int) # Start time cast into units of sample numbers
    acc_start_time = t_start*t_samp
    t_stop = np.ceil(stop_time/t_samp).astype(int) # Stop time cast into units of sample numbers
    startfile_index = np.where(cumsamples_per_file>=t_start)[0][0] # Index of first file that needs to be read
    stopfile_index = np.where(cumsamples_per_file>=t_stop)[0][0] # Index of last file that needs to be read
    startarray_index = ntsamples_per_file[startfile_index] + t_start - cumsamples_per_file[startfile_index] # Index of t_start in first data file to be read
    stoparray_index = ntsamples_per_file[startfile_index] + t_stop - cumsamples_per_file[stopfile_index] # Index of t_stop in last data file to be read
    print('No. of files to be read for given time range = %d'% (stopfile_index-startfile_index+1))

    for idx_file in range(startfile_index, stopfile_index+1):
        file_name = file_list[idx_file]
        if '/' in file_name:
            file_name = file_name.split('/')[-1]
        print('Reading in data from %s'% (file_name))
        data_file = fits.open(file_list[idx_file], ignore_missing_end=True)[-1].data['DATA'][:,:,:,:,0]
        # Data shape = (Nrows, NSBLK, NPOL, NCHANS)
        data_file = data_file.reshape((data_file.shape[0]*data_file.shape[1],data_file.shape[2],data_file.shape[3]))
        # Data shape = (Ntsamples, NPOL, NCHANS)
        data_file = np.moveaxis(data_file,0,-1) # Move time samples to last axis of array.
        # Select required polarization.
        data_file = data_file[pol] # Data shape = (len(pol), NCHANS, Ntsamples)
        # If only 1 file to load, ...
        if startfile_index==stopfile_index:
            return data_file[...,startarray_index:stoparray_index], hdr, acc_start_time
        # Keep track of data from multiple files.
        if idx_file==startfile_index:
            data = np.copy(data_file)[...,startarray_index:]
        elif idx_file==stopfile_index:
            data = np.concatenate((data, data_file[...,:stoparray_index]),axis=-1)
        else:
            data = np.concatenate((data, data_file),axis=-1)

    return data, hdr, acc_start_time
#################################################
