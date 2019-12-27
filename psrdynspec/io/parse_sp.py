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
    DMs = np.array(df['DM']) # DMs of single pulse candidates
    sigma = np.array(df['Sigma']) # Significance of single pulse candidate detections
    times = np.array(df['Time']) # Times at which single pulse candidates are found in the dedispersed time series
    samples = np.array(df['Sample']) # Sample numbers of candidates in dedispersed time series
    return DMs,sigma,times,samples
#########################################################################
# Apply a S/N cutoff to extracted candidates.
'''
Inputs:
DMs = Candidate DM (pc/cc)
sigma = Detection significance of candidate
times = Candidate times (s) in dedispersed time series
samples = Sample number of candidate in dedispersed time series.
sigma_min = S/N threshold
'''
def apply_sigma_cutoff(DMs,sigma,times,samples,sigma_min):
    indices = np.where(sigma>=sigma_min)[0]
    sigma_select = sigma[indices]
    DMs_select = DMs[indices]
    times_select = times[indices]
    samples_select = samples[indices]
    return DMs_select,sigma_select,times_select,samples_select
#########################################################################
# Remove duplicates from extracted candidate list.
'''
Inputs:
cand_DMs = DMs (pc/cc) of canidates
cand_sigma = Significance of candidate detection at respective DMs
cand_dedisp_times = Times (s) reported for candidates in dedispersed time dedisperesed time series
cand_dedisp_samples = Time sample numbers corresponding to above reported times
time_margin = A candidate identified within this margin (s) of another candidate at slightly different DM is marked as a duplicate detection.
'''
def remove_duplicates(cand_DMs,cand_sigma,cand_dedisp_times,cand_dedisp_samples,time_margin):
    # Sort the extracted candidate list in order of increasing times.
    sort_index = np.argsort(cand_dedisp_times)
    cand_dedisp_times = cand_dedisp_times[sort_index]
    cand_dedisp_samples = cand_dedisp_samples[sort_index]
    cand_DMs = cand_DMs[sort_index]
    cand_sigma = cand_sigma[sort_index]
    # Create new arrays for storing information of unique candidates.
    unique_cand_DMs = np.array([])
    unique_cand_sigma = np.array([])
    unique_cand_dedisp_times = np.array([])
    unique_cand_dedisp_samples = np.array([])
    # Remove duplicate candidates.
    i = 0
    while(i<len(cand_dedisp_times)):
        time = cand_dedisp_times[i]
        indices = i+np.where(np.abs(cand_dedisp_times[i:]-time)<=time_margin)[0]
        select_cand_index = indices[np.argmax(cand_sigma[indices])]
        # Save information about the selected candidate, if not duplicate.
        if (cand_dedisp_times[select_cand_index] not in unique_cand_dedisp_times):
            unique_cand_dedisp_samples = np.append(unique_cand_dedisp_samples,cand_dedisp_samples[select_cand_index])
            unique_cand_dedisp_times = np.append(unique_cand_dedisp_times,cand_dedisp_times[select_cand_index])
            unique_cand_DMs = np.append(unique_cand_DMs,cand_DMs[select_cand_index])
            unique_cand_sigma = np.append(unique_cand_sigma,cand_sigma[select_cand_index])
        # Update the iteration value by the number of indices parsed.
        i+= np.size(indices)
    print('Duplicate candidates removed.')
    return unique_cand_DMs, unique_cand_sigma, unique_cand_dedisp_times, unique_cand_dedisp_samples
#########################################################################
# Compile list of single pulse candidates over a specified DM range. Info from .singlepulse files is used.
'''
Inputs:
low_DM = Lowest DM (pc/cc) of interest
high_DM = Largest DM (pc/cc) of interest. Extracted single pulse candidates are identified at DMs between low_DM and high_DM (including endpoints).
glob_singlepulse = Glob search string to parse .singlepulse files
SINGLEPULSE_DIR = Path to .singlepulse files.
'''
def gen_singlepulse(low_DM,high_DM,glob_singlepulse,SINGLEPULSE_DIR):
    file_list = sorted(glob.glob(SINGLEPULSE_DIR+glob_singlepulse))
    file_DMs = np.array([float(file_name.split('DM')[-1].split('.singlepulse')[0]) for file_name in file_list])

    print('Collating details of candidate single pulses...')
    # Check which files correspond to DMs in the specified range.
    indices = np.where((file_DMs>=low_DM) & (file_DMs<=high_DM))[0]
    select_files = list(np.array(file_list)[indices])

    cand_DMs = np.array([]) # DM reported for candidate
    cand_sigma = np.array([]) # Significance of candidate detection.
    cand_dedisp_times = np.array([]) # Time (s) at which candidate is seen in the dedispersed time series
    cand_dedisp_samples = np.array([]) # Sample number of candidate in dedispersed time series

    # Generate full list of single pulse candidates from selected .singlepulse files
    for file_name in select_files:
        DMs,sigma,times,samples = read_singlepulse(file_name)
        cand_DMs = np.append(cand_DMs,DMs)
        cand_sigma = np.append(cand_sigma,sigma)
        cand_dedisp_times = np.append(cand_dedisp_times,times)
        cand_dedisp_samples = np.append(cand_dedisp_samples,samples)
    print('Collation complete.')
    return cand_DMs,cand_sigma,cand_dedisp_times,cand_dedisp_samples
#########################################################################
