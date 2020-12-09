# Parse .singlepulse files for single pulse candidates identified within a certain user-specified range in DM.
import numpy as np
import pandas as pd
import glob
#########################################################################
# Read in a .singlepulse file and organize its content into arrays.
'''
Input:
file_name = Name of input .singlepulse file to read (include path to file if not in current directory)
'''
def read_singlepulse(file_name):
    column_names = ['DM','Sigma','Time','Sample','Downfact'] # Column names in a .singlepulse file
    df = pd.read_csv(file_name,sep=r"\s+",header=None,names=column_names,comment='#')
    DMs = np.array(df['DM'],dtype='float64') # DMs of single pulse candidates
    sigma = np.array(df['Sigma'],dtype='float64') # Significance of single pulse candidate detections
    times = np.array(df['Time'],dtype='float64') # Times at which single pulse candidates are found in the dedispersed time series
    samples = np.array(df['Sample'],dtype='int') # Sample numbers of candidates in dedispersed time series
    downfact = np.array(df['Downfact'],dtype='int') # Smoothing factor to be applied on dedispersed time series for candidate detection
    return DMs, sigma, times, samples, downfact
#########################################################################
# Apply a S/N cutoff to extracted candidates.
'''
Inputs:
DMs = Candidate DM (pc/cc)
sigma = Detection significance of candidate
times = Candidate times (s) in dedispersed time series
samples = Sample number of candidate in dedispersed time series
downfact = Smoothing factor to be applied on dedispersed time series for candidate detection
sigma_min = S/N threshold
'''
def apply_sigma_cutoff(DMs,sigma,times,samples,downfact,sigma_min):
    indices = np.where(sigma>=sigma_min)[0]
    sigma_select = sigma[indices]
    DMs_select = DMs[indices]
    times_select = times[indices]
    samples_select = samples[indices]
    downfact_select = downfact[indices]
    return DMs_select,sigma_select,times_select,samples_select,downfact_select
#########################################################################
# Remove duplicates from extracted candidate list.
'''
Inputs:
cand_DMs = DMs (pc/cc) of canidates
cand_sigma = Significance of candidate detection at respective DMs
cand_dedisp_times = Times (s) reported for candidates in dedispersed time dedisperesed time series
cand_dedisp_samples = Time sample numbers corresponding to above reported times
time_margin = A candidate identified within this margin (s) of another candidate at slightly different DM is marked as a duplicate detection.
DM_margin = A candidate within this margin (pc/cc) of another candidate at a nearby time is marked as a duplicate detection.

Returns:
cand_DMs, cand_sigma, cand_dedisp_times, cand_dedisp_samples: Input arrays re-arranged in order of increasing time
select_indices = Indices of candidates selected for retention in output arrays. Indices correspond to those in output arrays sorted in ascending order of time.
'''
def remove_duplicates(cand_DMs,cand_sigma,cand_dedisp_times,cand_dedisp_samples,cand_downfact,time_margin,DM_margin):
    # Sort the extracted candidate list in order of increasing times.
    sort_index = np.argsort(cand_dedisp_times)
    cand_dedisp_times = cand_dedisp_times[sort_index]
    cand_dedisp_samples = cand_dedisp_samples[sort_index]
    cand_DMs = cand_DMs[sort_index]
    cand_sigma = cand_sigma[sort_index]
    cand_downfact = cand_downfact[sort_index]
    # Initialize empty arrays.
    check_indices = np.array([],dtype=int)
    select_indices = np.array([],dtype=int)
    # Remove duplicate candidates
    i = 0
    while i<len(cand_dedisp_times):
        if i not in check_indices:
            indices = np.where(np.logical_and(np.abs(cand_dedisp_times-cand_dedisp_times[i])<=time_margin, np.abs(cand_DMs-cand_DMs[i])<=DM_margin ))[0]
            select_cand_index = indices[np.argmax(cand_sigma[indices])]
            if select_cand_index not in check_indices:
                select_indices = np.append(select_indices, select_cand_index)
                check_indices = np.append(check_indices, indices)
        i+=1
    print('Duplicate candidates removed.')
    return cand_DMs, cand_sigma, cand_dedisp_times, cand_dedisp_samples, cand_downfact, select_indices
#########################################################################
# Compile list of single pulse candidates over a specified DM range. Info from .singlepulse files is used.
'''
Inputs:
low_DM = Lowest DM (pc/cc) of interest
high_DM = Largest DM (pc/cc) of interest. Extracted single pulse candidates are identified at DMs between low_DM and high_DM (including endpoints).
file_list = List of .singlepulse files to read
'''
def gen_singlepulse(low_DM,high_DM, file_list):
    file_DMs = np.array([float(file_name.split('DM')[-1].split('.singlepulse')[0]) for file_name in file_list])

    print('Collating details of candidate single pulses...')
    # Check which files correspond to DMs in the specified range.
    indices = np.where((file_DMs>=low_DM) & (file_DMs<=high_DM))[0]
    select_files = list(np.array(file_list)[indices])

    cand_DMs = np.array([]) # DM reported for candidate
    cand_sigma = np.array([]) # Significance of candidate detection.
    cand_dedisp_times = np.array([]) # Time (s) at which candidate is seen in the dedispersed time series
    cand_dedisp_samples = np.array([]) # Sample number of candidate in dedispersed time series
    cand_downfact = np.array([]) # Smoothing factor to be applied on dedispersed time series for candidate detection

    # Generate full list of single pulse candidates from selected .singlepulse files
    for file_name in select_files:
        DMs,sigma,times,samples,downfact = read_singlepulse(file_name)
        cand_DMs = np.append(cand_DMs,DMs)
        cand_sigma = np.append(cand_sigma,sigma)
        cand_dedisp_times = np.append(cand_dedisp_times,times)
        cand_dedisp_samples = np.append(cand_dedisp_samples,samples)
        cand_downfact = np.append(cand_downfact,downfact)
    print('Collation complete.')
    return cand_DMs,cand_sigma,cand_dedisp_times,cand_dedisp_samples,cand_downfact
#########################################################################
