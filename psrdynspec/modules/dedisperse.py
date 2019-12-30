# Set of functions for dedispersing dynamic spectra

# Import required packages.
import numpy as np
##########################################################################
# Calculate dispersive delay relative to a reference frequency for a given DM.
'''
Inputs:
freqs_GHz = 1D array of radio frequencies in GHz
DM = Dispersion measure (pc/cc) at which dispersive delay must be calculated
ref_freq = Reference frequency (GHz) relative to which dispersive delay must be calculated
'''
def calc_tDM(freqs_GHz,DM,ref_freq):
    kDM = 4.149e-3
    freq_factor = freqs_GHz**-2. - ref_freq**-2.
    tDM = kDM*DM*freq_factor
    return tDM
##########################################################################
# Brute-force dedisperse a dynamic spectrum at a given DM.
'''
Inputs:
ds = 2D array, dynamic spectrum
freqs_GHz = 1D array of frequencies (GHz) covered by dynamic spectrum
DM = DM at which dedispersion must be performed
ref_freq = Frequency (GHz) relative to which dispersive delay must be calculated
freq_low = Lowest frequency (GHz) for which we intend to calculate the dispersive delay
t_resol = Time resolution (s) of the data
start_time = Start time (s) of the data
'''
def dedisperse_ds(ds,freqs_GHz,DM,ref_freq,freq_low,t_resol,start_time):
    tDM = calc_tDM(freqs_GHz,DM,ref_freq)
    tDM_indices = np.round(tDM/t_resol).astype(int)

    # Determine no. of time slices to clip based on lowest frequency.
    max_tDM = calc_tDM(freq_low,DM,ref_freq)
    max_tDM_indices = np.round(max_tDM/t_resol).astype(int)

    n_chans,n_tsamples = ds.shape
    n_dedisp_tsamples = n_tsamples - max_tDM_indices
    dedisp_ds = np.zeros((n_chans,n_dedisp_tsamples))
    dedisp_times = start_time+np.arange(0,n_dedisp_tsamples)*t_resol
    for i in range(n_chans):
        start_t_index = tDM_indices[i]
        stop_t_index = n_dedisp_tsamples + tDM_indices[i]
        dedisp_ds[i] = ds[i,start_t_index:stop_t_index]
    dedisp_timeseries = np.nansum(dedisp_ds,axis=0)
    return dedisp_ds,dedisp_times,dedisp_timeseries
##########################################################################
# Dedisperse at a large number of trial DMs and determine the DM that maximizes S/N.
# Off-pulse RMS estimated using samples that lie outside of t_exclude seconds off the time series maximum.
'''
Inputs:
ds = 2D array, dynamic spectrum
freqs_GHz = 1D array of frequencies (GHz) covered by dynamic spectrum
trial_DMs = 1D array of trial DMs at which dedispersion must be performed
ref_freq = Frequency (GHz) relative to which dispersive delay must be calculated
freq_low = Lowest frequency (GHz) for which we intend to calculate the dispersive delay
t_resol = Time resolution (s) of the data
start_time = Start time (s) of the data
t_center = Central time (s) reported for candidate by PRESTO
t_exclude = No. of seconds of data around t_center should be exclude for off-pulse RMS calculation
'''
def calc_DM_at_maxSNR(ds,freqs_GHz,trial_DMs,ref_freq,freq_low,t_resol,start_time,t_center,t_exclude):
    # Keep track of how signal and offpulse_std vary with trial DM.
    signal_array = np.zeros_like(trial_DMs)
    offpulse_std_array = np.zeros_like(trial_DMs)
    # Determine which bins to ignore when computing off-pulse standard deviation of time series.
    t_exclude_bins = int(t_exclude//t_resol) # No. of bins to ignore on either side of pulse maximum
    for i in range(len(trial_DMs)):
        DM = trial_DMs[i]
        if (DM%2==0):
            print('Dedispersing at DM = %s'% (DM))
        dedisp_ds, dedisp_times,dedisp_timeseries = dedisperse_ds(ds,freqs_GHz,DM,ref_freq,freq_low,t_resol,start_time)
        signal = np.max(dedisp_timeseries) # Signal
        max_index = np.argmax(dedisp_timeseries) # Index corresponding to pulse maximum
        low_exclude_limit = max_index - t_exclude_bins # Lowest index to ignore
        high_exclude_limit = max_index + t_exclude_bins+1 # 1 greater than the highest index to ignore
        offpulse_std = np.std(np.append(dedisp_timeseries[:low_exclude_limit],dedisp_timeseries[high_exclude_limit:])) # Off-pulse standard deviation
        signal_array[i] = signal
        offpulse_std_array[i] = offpulse_std
    # Find the DM at which S/N gets maximized.
    index_max_signal = np.argmax(signal_array)
    offpulse_std_value = offpulse_std_array[index_max_signal]
    SNR = signal_array/offpulse_std_value
    return SNR, signal_array, offpulse_std_value
##########################################################################
# Fuction to calculate block length (samples) given total number of time samples (T), no. of blocks (N) and dispersive time samples (t_DM)
'''
Inputs:
tot_time_samples = No. of time samples in time series
N_blocks = No. of blocks in which data are to be loaded
max_tDM_samples = Max. no. of time samples covered by dispersive sweep across the observing band
'''
def calc_block_length(tot_time_samples,N_blocks,max_tDM_samples):
    block_length = np.ceil((tot_time_samples+(N_blocks-1)*max_tDM_samples)/N_blocks).astype(int)
    return block_length
##########################################################################
