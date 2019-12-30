# Load data from list of PSRFITS files.

from astropy.io import fits
import glob
import numpy as np
#################################################
# Extract data from PSRFITS file.
'''
Inputs:
glob_fits_files = Glob string to parse .fits files
FITS_DIR = Path to .fits files
pol = List of polarization indices to extract (default = [0])
'''
def load_psrfits_data(glob_fits_files,FITS_DIR,pol=[0]):
    file_list = sorted(glob.glob(FITS_DIR+glob_fits_files))
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
