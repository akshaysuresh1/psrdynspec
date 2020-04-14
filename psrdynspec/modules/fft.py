# Tools for performing FFT of a signal

import numpy as np
from scipy.signal import find_peaks
#############################################################
# Perform |FFT|^2 of 1D data blocks of length Nfft and average them to obtain power spectrum.
# For a single N-point transform, set Nfft = np.size(data_1d).
'''
Inputs:
data1d = 1D array of data values
Nfft = Length of FFT for each block
tsamp  = Time sampling (s) of data
remove_DC_spike = True/False
return_positive_freqs_only = True/False
'''
def fft1d_blocks_avg(data1d, Nfft, tsamp, remove_DC_spike=False,return_positive_freqs_only=False):
    T = len(data1d) # No. of time samples in data
    Nfft = int(Nfft) # Cast Nfft as integer ignoring decimal part.
    if Nfft%2!=0:
        print('Input FFT length is not divisible by 2. Rounding off to nearest multiple of 2 greater than input length.')
        Nfft = int(np.ceil(Nfft/2)*2)
    remainder = T % Nfft
    if remainder!=0:
        print('Padding data with median values to the nearest size divisibly by %d'% (Nfft))
        pad_size = int(Nfft - remainder)
        data1d = np.concatenate([data1d, np.ones(pad_size)*np.nanmedian(data1d)])
    N_windows = np.ceil(2*T/Nfft).astype(int) - 1 # No. of STFT windows. Adjacent windows are offset by half a window length.
    print('Total number of blocks to FFT = ',N_windows)
    print('FFT length = ',Nfft)
    half_window_length = Nfft//2
    # Initialize array for storing power spectrum values.
    power_spectrum = np.zeros(Nfft)
    # Perform block-wise FFT and sum power spectra over blocks.
    for n in range(N_windows):
        if (np.mod(n+1,N_windows/4)==0): print("Processing block", n+1)
        xb = data1d[n*half_window_length:n*half_window_length+Nfft]
        if remove_DC_spike:
            xb -= np.mean(xb)
        xfft = np.fft.fft(xb)
        spec = np.abs(xfft)**2/Nfft
        power_spectrum += spec
    power_spectrum = power_spectrum/N_windows
    # Reorder frequencies. Bring the zero frequency to the center of the array.
    power_spectrum = np.fft.fftshift(power_spectrum)
    # Spectral resolution achieved with Nfft-length FFT
    freq_resol = 1./(Nfft*tsamp)
    freq_array = np.arange(-Nfft/2,Nfft/2)*freq_resol   # Hz
    # Return only positive frequencies if specified.
    if return_positive_freqs_only:
        pos_indices = np.where(freq_array>0)[0]
        freq_array = freq_array[pos_indices]
        power_spectrum = power_spectrum[pos_indices]
    return power_spectrum, freq_array
#############################################################
# Search for peaks in a power spectrum.
'''
Inputs:
frequencies = 1D array of Fourier frequency values (Hz)
power_spectrum = 1D array of power_spectrum values
search_scale = 'linear' or 'log' (d: log)
niter = No. of iterations (integer) to perform for iterative sigma clipping (d: 5)
sigma_clip = No. of standard deviations outside of median to be clipped (d: 3.0)
'''
def find_FFTpeaks(frequencies,power_spectrum,search_scale='log',niter=5,sigma_clip=3.0):
    if search_scale=='log':
        mask_array = np.log(power_spectrum)
    elif search_scale=='linear':
        mask_array = power_spectrum
    else:
        return "Search scale not defined. Choose either 'linear' or 'log'."
    median = np.median(mask_array)
    std = np.std(mask_array)
    print('Iteratively masking elements within 3 sigma of median in %s scale'% (search_scale))
    for i in range(niter):
      mask_array = np.ma.masked_outside(mask_array,median-sigma_clip*std,median+sigma_clip*std)
      median = np.median(mask_array)
      std_new = np.std(mask_array)
      std_diff_percent = np.abs(std - std_new)*100/std_new # Percentage change in standard deviation between iterations.
      if (std_diff_percent<=5):
          print('Convergence criterion achieved at iteration %d'% (i+1))
          break
      elif (i==niter-1):
          print('Max iterations limit reached.')
    if search_scale=='log':
        median = np.exp(median)
        std = np.exp(std)
    peak_indices, _ = find_peaks(power_spectrum,height=median+5*std)
    peak_freqs = frequencies[peak_indices]
    peak_powerspec = power_spectrum[peak_indices]
    return peak_indices, peak_freqs, peak_powerspec
#############################################################
# Peform FFT with noise masking, if specified.
'''
Inputs:
data1d = 1D array of data values
Nfft = FFT length (preferably a power of 2)
axis_resol = Resolution along axis element of data aray. For example, time resolution if data_1d is a time series
N_sigma = Factor N to be used for sigma clipping. If exclude_noise=True, values less than N*sigma of the median are replaced with zeros before FFT computation.
exclude_noise = Do you want to exclude the baseline statistical noise when performing FFT? (True/False) (d: 3.0)
search_scale = 'linear' or 'log' (d: log)
niter = No. of iterations (integer) to perform for iterative sigma clipping (d: 5)
sigma_clip = No. of standard deviations outside of median to be clipped (d: 3.0)
'''
def fft1d_mask(data1d,Nfft,axis_resol,exclude_noise=False,N_sigma=3.0,search_scale='log',niter=5,sigma_clip=3.0):
    # If specified, mask noise before computing FFT.
    if (exclude_noise==True):
        threshold = np.median(data1d)+N_sigma*np.std(data1d)
        data1d = np.ma.filled(np.ma.masked_less(data1d ,threshold),0) # Replace values below the threshold with zeros.
    # Compute FFT in blocks on length Nfft.
    power_spectrum, frequencies = fft1d_blocks_avg(data1d, Nfft, axis_resol, remove_DC_spike=True,return_positive_freqs_only=True)
    peak_indices, peak_freqs, peak_powerspec = find_FFTpeaks(frequencies,power_spectrum,search_scale,niter,sigma_clip)
    return frequencies, power_spectrum, peak_indices, peak_freqs, peak_powerspec
#############################################################
