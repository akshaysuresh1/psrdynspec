# Load data from list of PSRFITS files.

from astropy.io import fits
import numpy as np
import sys
#################################################
# Extract data from PSRFITS file.
'''
Inputs:
file_list = List of .fits files
FITS_DIR = Path to .fits files
pol = List of polarization indices to extract (default = [0])
'''
def load_psrfits_data(file_list,FITS_DIR,pol=[0]):
    if len(file_list)==0:
        print('PSRFITS file list is empty. Quitting.')
        sys.exit(1)
    for i in range(len(file_list)):
        file_name = file_list[i].split('/')[-1]
        print('Reading in data from %s'% (file_name))
        data_file = fits.open(file_list[i])[-1].data['DATA'][:,:,:,:,0]
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
