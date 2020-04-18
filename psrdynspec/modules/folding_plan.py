'''
NOTE: This script is an adapted version of https://github.com/v-morello/riptide/blob/master/riptide/processing_plan.py.

MIT License

Copyright (c) 2017-2020 Vincent Morello

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''
import numpy as np
import pandas as pd
#################################################################################
# FUNCTION DEFINITIONS
'''
# Return a pandas DataFrame with the processing plan for folding searches of a timeseries.
# Period ranges get grouped in octaves of form [2^k Pmin, 2^(k+1) * Pmin).

Inputs:
P_min = Minimum period (s) to search
P_max = Maximum period (s) to search
bins_min = Minimum no. of bins, Pmin/bins_min sets the minimum pulse width resolution.
tau = Sampling time (s) in input timeseries

Returns:
octaves = pandas DataFrame object containing processing details of each octave
'''
def gen_octaves(P_min, P_max, bins_min, tau):
    octaves = []
    column_names = ['dsfactor','tsamp','bins_min','bins_max','period_min','period_max']

    # Cast input to correct types.
    bins_min = int(bins_min)
    tau = float(tau)

    if tau * bins_min > P_min:
        raise ValueError('Pmin must be larger than tau x bins_min.')

    # Downsampling factor for first octave.
    ds = int(np.floor(P_min/(bins_min*tau)))
    bins_max = 2*bins_min # Increase no. of bins by a factor of 2 from start to end of octave period range.

    while True:
        tsamp = ds*tau
        min_octave_period = bins_min*tsamp
        max_octave_period = 2*min_octave_period
        current_octave = [ds, tsamp, bins_min, bins_max, min_octave_period, max_octave_period]
        octaves.append(current_octave)

        # Downsample further by a factor of 2 on moving from one octave to the next.
        # This keeps bins_min preserved across octaves.
        ds *= 2

        if (max_octave_period>=P_max):
            break

    octaves = pd.DataFrame(octaves, columns=column_names)

    # Edit last octave to stop folding search at a period sufficiently close to Pmax.
    last_octave = octaves.iloc[-1]
    last_bins_max = int(np.ceil(P_max/last_octave.tsamp))
    last_periods_max = last_bins_max*last_octave.tsamp

    index = len(octaves)-1
    octaves.loc[index,'bins_max'] = last_bins_max
    octaves.loc[index,'period_max'] = last_periods_max

    return octaves

'''
Generate 1D array of trial periods based on an input processing plan from gen_octaves(...).

Inputs:
octaves = pandas DataFrame object containing processing details of each octave
tau = Sampling time (s) in input timeseries to fold at multiple periods
N = No. of samples in input timeseries

Returns:
periods = 1D array of trial folding periods (s)
fold_bins = 1D array of phase bins to use for folding at above trial periods
'''
def gen_periods(octaves, tau, N):
    periods = [] # Construct 1D array of periods assuming optimal spacing from FFA search.
    fold_bins = [] # 1D array of phase bins to use for folding at above trial periods.
    for _, step in octaves.iterrows():
        ds = step.dsfactor
        ns = int(N/ds)
        bins = np.arange(int(step.bins_min), int(step.bins_max))
        for b in bins:
            m = int(ns/b)
            fold_bins.append(np.ones(m)*b)
            for sh in range(m):
                periods.append(ds*b + ds*sh*b/(m*b - sh))
    periods = np.array(periods)*tau
    fold_bins = np.hstack(np.array(fold_bins))

    # A period at the final slope s=(m-1) of a base period p (samples) will be covered at a small value of s for the next base period (p+1)..
    # Remove instances of such redundant periods.
    mask = np.ones(periods.size, dtype=bool)
    # Reverse 1D array of periods and mask values such that the masked array is strictly decreasing.
    rev_periods = periods[::-1]
    diff = np.diff(rev_periods)

    # Find indices where the immediately succesive element is greater.
    indices = np.where(diff > 0)[0]
    for ix in indices:
        for jj in range(ix+1, rev_periods.size):
            mask[jj] = False
            if rev_periods[jj] < rev_periods[ix]:
                break
    mask = mask[::-1] # Reverse mask.            

    # Find indices where the trial periods are greater than the maximum period of the last octave.
    indices = np.where(periods>=np.max(np.array(octaves.period_max)))[0]
    if len(indices)!=0:
        mask[indices] = False

    return periods[mask], fold_bins[mask]

'''
Updates pandas DataFrame to include no. of periods covered in each octave.

Inputs:
octaves = pandas DataFrame object containing processing details listed octave-wise,
periods = 1D array of trial folding periods (s)
'''
def update_processing_plan(octaves, periods):
    column_names = octaves.columns.tolist()
    column_names.append('N_periods')

    steps_updated = []
    for _, step in octaves.iterrows():
        # Perform step-wise retrieveal of values of different quantities in processing plan.
        ds = int(step.dsfactor)
        tsamp = float(step.tsamp)
        bins_min = int(step.bins_min)
        bins_max = int(step.bins_max)
        period_min = float(step.period_min)
        period_max = float(step.period_max)
        N_periods_in_octave = len(np.where( np.logical_and(periods>=period_min, periods<period_max))[0])
        # Updated step properties
        current_step = [ds, tsamp, bins_min, bins_max, period_min, period_max, N_periods_in_octave]
        steps_updated.append(current_step)

    steps_updated = pd.DataFrame(steps_updated, columns=column_names)

    return steps_updated
#################################################################################
# CLASS DEFINITIONS
class ProcessingPlan(object):
    ''' Defines the sequence of downsamplings and octave searches to perform on a time series.
    '''

    def __init__(self, nsamp, tsamp, bins_min, octaves, periods, fold_bins):
        '''
        This should not be called directly. Use ProcessingPlan.create() instead.
        '''
        self.nsamp = nsamp # No. of samples in input timeseries
        self.tsamp = tsamp # Sampling time (s) of input timeseries
        self.bins_min = bins_min # Minimum no. of bins
        self.periods = periods
        self.octaves = octaves
        self.period_min = np.min(periods)
        self.period_max = np.max(periods)
        self.fold_bins = fold_bins
        self.min_width_res = self.period_min/self.bins_min # Minimum pulse width resolution

    @staticmethod
    def create(nsamp, tsamp, bins_min, P_min, P_max):
        octaves = gen_octaves(P_min, P_max, bins_min, tsamp)
        periods, fold_bins = gen_periods(octaves, tsamp, nsamp)
        octaves = update_processing_plan(octaves, periods)
        return ProcessingPlan(nsamp, tsamp, bins_min, octaves, periods, fold_bins)

    def __repr__(self):
        lines = ['Properties of input timeseries:',
            'Number of samples  = %d' % (self.nsamp),
            'Sampling time (s)  = %.6e' % (self.tsamp),
            'Min. no. of bins   = %d'% (self.bins_min),
            'Min. pulse width resolution (s) = %.6e'% (self.min_width_res),
            '----------------------------------------------------------------',
            '',
            repr(self.octaves)
            ]
        return '\n'.join(lines)
#################################################################################
