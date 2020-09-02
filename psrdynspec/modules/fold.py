# Functions to fold a timeseries or dynamic spectrum at given period, or generate a processing plan for bulk folding.
# Only brute force time-domain folding supported as of now.
import numpy as np
import pandas as pd
from psrdynspec.modules.filters1d import blockavg1d
from astropy.stats import sigma_clip
##########################################################################
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
    phibin_edges = np.linspace(0.,1.,Nbins+1) # Edges of phase bins
    phibin_centers = 0.5*(phibin_edges[1:]+phibin_edges[:-1]) # Centers of phase bins
    indbin_edges = np.digitize(phi,phibin_edges,right=True) # Place phase values in phase bins. Returns indices of placement in phase bins array.
    # Given edges x < y, a phase bin is defined to be the range (x,y].

    # Arrays to store profile and counts per phase bin.
    profile = np.zeros(Nbins)
    counts = np.zeros(Nbins)
    for n in range(1,Nbins+1):
        profile[n] = np.sum(timeseries[np.where(indbin_edges==n)]) # Sum up time series flux values that fall in one phase bin.
        counts[n] = np.size(np.where(indbin_edges==n)) # No. of counts in bin.
    profile /= counts # Divide by the number of counts to generate an average profile.
    return profile, phibin_centers
##########################################################################
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
##########################################################################
# Select the metric function corresponding to specified metric string.
'''
Inputs:
metric = Metric to maximize for determining optimal folding period ('reducedchisquare', 'profmax', 'profSNR')
'''
def select_metric_function(metric):
    if (metric=='profmax'):
        metric_function = calc_profilemax
    elif (metric=='reducedchisquare'):
        metric_function = calc_reduced_chisquare_profile
    elif (metric=='profSNR'):
        metric_function = calc_sn
    else:
        raise Exception('Chosen metric not recognized.')
        return None
    return metric_function
##########################################################################
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
    metric_function = select_metric_function(metric)
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
##########################################################################
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
def calc_reduced_chisquare_profile(profile):
    profile_without_outliers = sigma_clip(profile, sigma=5.0, cenfunc='median', stdfunc='std', maxiters=3)
    std_offpulse = np.std(profile_without_outliers)
    med = np.ma.median(profile_without_outliers)
    reduced_chisquare = np.mean((profile-med)**2)/(std_offpulse**2)
    return reduced_chisquare
###########################################################################
# Calculate signal-to-noise ratio of a profile.
'''
Inputs:
profile = Folded pulse profile
'''
def calc_sn(profile):
    profile_without_outliers = sigma_clip(profile, sigma=3.0, cenfunc='median', stdfunc='std', maxiters=5)
    noise = np.std(profile_without_outliers)
    signal = np.max(profile)
    sn = signal/noise
    return sn
###########################################################################
# Execute a folding search on a timeseries based on a plan contained in a ProcessingPlan(..) object.
'''
Inputs:
timeseries = 1D timeseries to fold
times = 1D array of times (seconds)
plan = a ProcessingPlan(..) object describing sequence of downsampling and partial searches to be performed
metric = 'reduced chi square', 'profmax' : Metric for determining optimal folding period.
'''
def execute_plan(timeseries, times, plan, metric):
    metric_function = select_metric_function(metric)
    periods = plan.periods # 1D array of trial periods
    fold_bins = plan.fold_bins.astype(int) # 1D array of phase bins to use for folding at above trial periods
    metric_values = np.zeros(len(periods)) # 1D array to store metric values for above trial periods
    print('Folding data at a large number of trial periods...')
    print('Metric: ',metric)
    print('Index      Period (s)       Metric value')
    count = 0 # Keep track of cumulative number of periods covered after each octave.
    dsfactor_values = np.zeros(len(periods)) # Keep track of downsampling factors for different trial periods.
    for _, step in plan.octaves.iterrows():
        dsfactor = int(step.dsfactor)  # Downsampling factor
        N_periods = int(step.N_periods) # No. of periods covered in each octave
        period_min = float(step.period_min) # Minimum octave trial period (s)
        period_max = float(step.period_max) # Maximum octave trial period (s)
        print('Downsampling input time series by factor %d for octave periods %.4f - %.4f s'% (dsfactor, period_min, period_max))
        blkavg_timeseries = blockavg1d(timeseries, dsfactor) # Downsample timeseries by above factor.
        blkavg_times = blockavg1d(times, dsfactor) # Block average 1D array of times by same factor.
        for i in range(N_periods):
            index = count+i
            profile, phibins = fold_ts(blkavg_timeseries, blkavg_times, periods[index], int(fold_bins[index]))
            metric_values[index] = metric_function(profile)
            print('%3d         %8.6f        %8.3f'% (index, periods[index], metric_values[index]))
        dsfactor_values[count:count+N_periods] =  dsfactor
        count += N_periods

    # Find trial period that maximizes the chosen metric.
    global_metricmax_index = np.nanargmax(metric_values)
    global_metricmax = np.nanmax(metric_values)
    best_period = periods[global_metricmax_index]
    optimal_bins = fold_bins[global_metricmax_index]
    optimal_dsfactor = dsfactor_values[global_metricmax_index]
    print('Metric is maximized is at index %d, i.e., P = %8.6f s.'% (global_metricmax_index,best_period))

    return metric_values, global_metricmax_index, global_metricmax, best_period, optimal_bins, optimal_dsfactor
###########################################################################
