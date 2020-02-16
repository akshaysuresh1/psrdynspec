# Read in a PRESTO dedispersed time series (.dat file)
import numpy as np
##############################################################
# Extract dedispersed time series from a PRESTO .dat file.
'''
Inputs:
dat_file = Name of .dat file to load
'''
def load_presto_dat(dat_file):
    dedispersed_timeseries = np.fromfile(dat_file,dtype=np.float32)
    return dedispersed_timeseries
##############################################################
