import numpy as np

############################################################################
# Construct the full dynamic spectrum from a collection of subbands.
# Gaps in frequency coverage are filled with NaNs.
'''
Inputs:
data_shape = Shape of raw dynamic spectrum
subband_data = List of 2D arrays, len(list) = No. of subbands
freqs_subband = List of 1D arrays, 1D arrays specify frequencies covered in respective subbands
freqs_GHz = 1D array of radio frequencies (GHz)
'''
def construct_ds_from_subbands(data_shape,subband_data,freqs_subband,freqs_GHz):
    whole_ds = np.zeros(data_shape)+np.NaN
    # Keep track of the lowest and highest frequencies of the DS that are not flagged.
    min_freq_index = np.size(whole_ds)-1
    max_freq_index = 0
    for i in range(len(freqs_subband)):
        in1 = np.where(freqs_GHz==freqs_subband[i][0])[0][0]
        in2 = np.where(freqs_GHz==freqs_subband[i][-1])[0][0]
        min_freq_index = np.min([min_freq_index,in1])
        max_freq_index = np.max([max_freq_index,in2])
        whole_ds[in1:in2+1] = subband_data[i]
    whole_ds = whole_ds[min_freq_index:max_freq_index+1]
    freqs = freqs_GHz[min_freq_index:max_freq_index+1]
    return whole_ds, freqs
############################################################################
