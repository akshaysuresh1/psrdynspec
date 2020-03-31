# Tools for performing FFT of a signal

import numpy as np
from scipy.signal import find_peaks
#############################################################
# Perform |FFT|^2 of 1D data blocks of length Nfft and sum them to obtain power spectrum.
# For a single N-point transform, set Nfft = np.size(data_1d).
'''
Inputs:
data_1d = 1D array of data values
Nfft = Length of FFT for each block
t_resol  = Time resolution (s) of data
remove_DC_spike = True/False
'''
def fft1d_blocks_avg(data_1d,Nfft,t_resol,remove_DC_spike=False):
    N_samples = np.size(data_1d)
    print('Total number of samples in the data = ',N_samples)
    if (N_samples<Nfft):
        print('Input array size < FFT length')
        print('Padding input array with median values to the specified FFT length')
        data_1d = np.append(data_1d, np.ones(Nfft-N_samples)*np.nanmedian(data_1d))
    # Total number of independent blocks to FFT
    Nblocks_max = N_samples//Nfft
    print('Total number of blocks to FFT = ',Nblocks_max)
    print('FFT length = ',Nfft)
    # Initialize array for storing power spectrum values.
    power_spectrum = np.zeros(Nfft)
    # Perform block-wise FFT and sum power spectra over blocks.
    for n in range(Nblocks_max):
        if (np.mod(n+1,Nblocks_max/4)==0): print("Processing block", n+1)
        xb = data_1d[n*Nfft:(n+1)*Nfft]
        if (remove_DC_spike==True):
            xb -= np.mean(xb)
        xfft = np.fft.fft(xb)
        spec = np.abs(xfft)**2/Nfft
        power_spectrum += spec
    # Reorder frequencies. Bring the zero frequency to the center of the array.
    power_spectrum = np.fft.fftshift(power_spectrum)
    # Spectral resolution achieved with Nfft-length FFT
    freq_resol = 1./(Nfft*t_resol)
    freq_array = np.arange(-Nfft/2,Nfft/2)*freq_resol
    return power_spectrum,freq_array
#############################################################
# Peform FFT with noise masking, if specified.
'''
Inputs:
data1d = 1D array of data values
Nfft = FFT length (preferably a power of 2)
axis_resol = Resolution along axis element of data aray. For example, time resolution if data_1d is a time series
N_sigma = Factor N to be used for sigma clipping. If exclude_noise=True, values less than N*sigma of the median are replaced with zeros before FFT computation.
exclude_noise = Do you want to exclude the baseline statistical noise when performing FFT? (True/False)
'''
def fft1d_mask(data1d,Nfft,axis_resol,N_sigma,exclude_noise=False):
    # If specified, mask noise before computing FFT.
    if (exclude_noise==True):
        threshold = np.median(data1d)+N_sigma*np.std(data1d)
        data1d = np.ma.filled(np.ma.masked_less(data1d ,threshold),0) # Replace values below the threshold with zeros.
    # Compute FFT in blocks on length Nfft.
    power_spectrum, frequencies = fft1d_blocks_avg(data1d,Nfft,axis_resol,remove_DC_spike=True)
    # Only consider the positive frequencies.
    positive_freqs = frequencies[np.where(frequencies>0)]
    powerspec_at_posfreqs = power_spectrum[np.where(frequencies>0)]
    # Calculate peaks in the power spectrum.
    peak_indices = find_peaks(powerspec_at_posfreqs)[0] # Indices of peaks in power spectrum.
    peak_freqs = positive_freqs[peak_indices] # Peak frequencies of power power spectrum
    peak_powerspec_values = powerspec_at_posfreqs[peak_indices] # Power spectrum values at peak frequencies
    return positive_freqs, powerspec_at_posfreqs, peak_indices, peak_freqs, peak_powerspec_values
#############################################################
