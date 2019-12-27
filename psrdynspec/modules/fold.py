# Functions to fold a timeseries or dynamic spectrum at given period.

import numpy as np
##############################################################
# Fold a time series at a given period.
'''
Inputs:
timeseries = 1D time series to fold
times = 1D array of times (s)
pfold = Folding period (s)
Nbins = No. of phase bins
'''
def fold_ts(timeseries, times, pfold, Nbins):
    phi = times/pfold % 1.0
    phibins = np.linspace(0.,1.,Nbins)
    indbins = np.digitize(phi,phibins,right=True) # Place phase values in phase bins. Returns indices of placement in phase bins array.

    # Arrays to store profile and counts per phase bin.
    profile = np.zeros(Nbins)
    counts = np.zeros(Nbins)
    for n in range(Nbins):
        profile[n] = np.sum(timeseries[np.where(indbins==n)]) # Sum up timeseries values that belong to one phase bin.
        counts[n] = np.size(np.where(indbins==n)) # No. of counts in bin.
    profile[1:] /= counts[1:] # Divide by the number of counts to generate an average profile.
    return profile[1:], phibins[1:]
##############################################################
# Convert a long time series into number of rotations and phase bins based on a given rotation period.
'''
Inputs:
timeseries = 1D time series to fold
times = 1D array of times (s)
pfold = Folding period (s)
Nbins = No. of phase bins
'''
def fold_rotations_ts(timeseries, times, pfold, Nbins):
    rotations = times//pfold
    N_rotations = int(np.max(rotations)+1) # Number of rotations including rotation 0.
    print('Total no. of rotations in data = %d'% (N_rotations))

    phi = times/pfold % 1.0
    phibins = np.linspace(0.,1.,Nbins)
    indbins = np.digitize(phi,phibins,right=True) # Place phase values in phase bins. Return indices of placement in phase bins array.

    profile_rotations = np.zeros((N_rotations,Nbins-1))
    counts_perrot_phibin = np.zeros((N_rotations,Nbins-1)) # Counts per rotation per phase bin
    for i in range(N_rotations):
        print ('Processing rotation %d'% (i))
        data_indices = np.where(rotations==i)[0] # Indices in the timeseries that belong to rotation i
        data_rotation = timeseries[data_indices] # Time series for rotation i
        phase_bin_indices = indbins[data_indices] # Phase bin indices for data elements in rotation i

        profile = np.zeros(Nbins)
        counts = np.zeros(Nbins)
        for n in range(Nbins):
            profile[n] = np.sum(data_rotation[np.where(phase_bin_indices==n)]) # Sum up timeseries values that belong to one phase bin.
            counts[n] = np.size(np.where(phase_bin_indices==n)) # No. of counts in bin.
        profile[1:] /= counts[1:] # Divide by the number of counts to generate an average profile.
        profile_rotations[i] = profile[1:]
        counts_perrot_phibin[i] = counts[1:]
    #profile_rotations -= np.nanmedian(profile_rotations.flatten())# Set baseline flux to zero.
    return profile_rotations, counts_perrot_phibin, phibins[1:]
##############################################################
# Fold at a large number of specified periods and keep track of the chosen metric for determining optimal folding period.
'''
Inputs:
timeseries = 1D time series to fold
times = 1D array of times (s)
periods = 1D array of trial folding period (s)
Nbins = No. of phase bins
metric = 'reduced chi square', 'profmax' : Metric for determining optimal folding period.
'''
def fold_metric_periods(timeseries,times,periods,Nbins,metric):
    metric_values = np.zeros(len(periods))
    if (metric=='profmax'):
        metric_function = calc_profilemax
    elif (metric=='reduced chi square'):
        metric_function = calc_reduced_chiquare_profile
    print('Folding data at a large number of trial periods...')
    print('Metric: ',metric)
    print('Index      Period (s)       Metric value')
    for i in range(len(periods)):
        pfold = periods[i]
        profile,phibins = fold_ts(timeseries, times, pfold, Nbins)
        metric_values[i] = metric_function(profile)
        print('%3d         %8.6f        %8.3f'% (i, pfold,metric_values[i]))

    global_metricmax_index = np.nanargmax(metric_values)
    global_metricmax = np.nanmax(metric_values)
    best_period = periods[global_metricmax_index]
    print('Metric is maximized is at index %d, i.e., P = %8.6f s.'% (global_metricmax_index,best_period))
    return metric_values, global_metricmax_index, global_metricmax, best_period
##############################################################
# Calculate profile maximum.
'''
Inputs:
profile = Folded pulse profile
'''
def calc_profilemax(profile):
    return np.max(profile)
##########################################################################
# Calculate reduced chi square of a pulse profile relative to noise.
'''
Inputs:
profile = Folded pulse profile
'''
def calc_reduced_chiquare_profile(profile):
    med = np.median(profile)
    profile_without_outliers = np.ma.masked_outside(profile,med-3*np.std(profile),med+3*np.std(profile))
    std_offpulse = np.std(profile_without_outliers)
    reduced_chisquare = np.mean((profile-med)**2)/(std_offpulse**2)
    return reduced_chisquare
###########################################################################
