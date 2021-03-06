# Module to read data from a filterbank file.
import numpy as np
##############################################################
# Extract data from filterbank file.
'''
Inputs:
f = File object with cursor initialized inside a filterbank file (e.g, f = open(fil_file, 'rb'))
t_start = Sample number corresponding to the start time
t_stop = Sample number corresponding to the stop time
n_ifs = No. of polarizations in .fil file
nchans = Total no. of spectral channels
n_bytes = No. of bytes per data sample
hdr_size = Header size (bytes)
pol = List of polarization indices to extract (default = [0])
current_cursor_position = Current position of the file cursor (default = 0)
'''
def load_fil_data(f,t_start,t_stop,n_ifs,nchans,n_bytes,hdr_size,pol=[0],current_cursor_position=0):
    datatype = setup_dtype(n_bytes)
    N_t_samples = int(t_stop-t_start)
    start_position = int(t_start*n_ifs*nchans*n_bytes)+hdr_size # Position where file cursor must be placed before reading.
    f.seek(start_position-current_cursor_position,1)
    data = np.fromfile(f,count=int(nchans*n_ifs*N_t_samples),dtype=datatype)
    data = data.reshape((N_t_samples,n_ifs,nchans))
    data = np.moveaxis(data,0,-1)[pol]
    return data
##############################################################
# Set up data type for data array according to the no. of bits per sample.
'''
Inputs:
n_bytes = No. of bytes per pixel of the dynamic spectrum.
'''
def setup_dtype(n_bytes):
    if n_bytes==8:
        return 'float64'
    elif n_bytes == 4:
        return 'float32'
    elif n_bytes == 2:
        return 'uint16'
    elif n_bytes == 1:
        return 'uint8'
    else:
        print('Having trouble setting dtype, assuming float32.')
        return 'float32'
##############################################################
