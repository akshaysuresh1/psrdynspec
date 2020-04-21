# Various normalizations applied to a 1D array (timeseries / power spectra)

import numpy as np
from astropy.stats import sigma_clip
#####################################################################
# Normalize 1D array to zero mean and unit variance.
'''
Inputs:
array1d = 1D array of values
'''
def normalize_stdnormal(array1d):
    clip_array = sigma_clip(array1d,sigma=3,cenfunc='median', stdfunc='std')
    mean = np.mean(clip_array)
    sigma = np.std(clip_array)
    array1d = (array1d - mean)/sigma
    return array1d
#####################################################################
# Normalize the maximum of a 1D array to unity.
'''
Inputs:
arrayNd = A Nd-array of values
'''
def normalize_unitpeak(arrayNd):
    normalization_factor = 1./np.nanmax(arrayNd)
    norm_array = normalization_factor*arrayNd
    return norm_array
#####################################################################
