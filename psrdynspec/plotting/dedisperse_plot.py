from .config import *
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec

# Plots saved to disk by default. Plots displayed on live window only if show_plot==True.
##########################################################################
# Plot S/N vs. DM for a single pulse candidate.
'''
Inputs:
trial_DMs = 1D array of trial DMs (pc/cc) at which dedispersion must be performed
SNR = 1D array of signal-to-noise ratios at above trial DMs
optimal_DM = Optimal DM according to prescribed metric
t_cand = Time (s) of occurrence of candidate
basename = Basename (string) for output plot (including path)
show_plot = Do you want to show the plot live? (True/False) (default = False)
'''
def plot_SNRvsDM(trial_DMs,SNR,optimal_DM,t_cand,basename,show_plot=False):
    plot_name = basename+'_t'+'%07.3f'% (t_cand)+'_SNRvsDM.png'
    print('Plotting S/N vs. trial DM for candidate at t = %.3f s'% (t_cand))
    plt.plot(trial_DMs,SNR)
    plt.xlabel('Trial DM (pc cm$^{-3}$)',fontsize=14)
    plt.ylabel('S/N',fontsize=14)
    plt.axvline(x=optimal_DM,label='DM$_{\mathrm{opt}}$',linestyle='--',color='r')
    plt.legend(loc='best',prop={'size':12})
    plt.savefig(plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
# In a single figure, combine plots of a dedispersed dynamic spectrum subband, a dedispersed time series and S/N vs. trial DM.
'''
Inputs:
dedisp_ds = 2D array, dedispersed dynamic spectrum
dedisp_timeseries = Sum of dedisp_ds over frequency
dedisp_times = 1D array of times(s) covered in dedisp_timeseries
freqs_array = 1D array of frequencies (GHz) covered in dedisp_ds
t_cand = Time (s) of occurrence of candidate
t_resol  = Time resolution (s) of data
trial_DMs = 1D array of trial DMs at which dedispersion must be performed
SNR = 1D array of signal-to-noise ratios at above trial DMs
optimal_DM = DM that maximized S/N
offpulse_std_value = Off-pulse standard deviation to normalize the dedispersed time series.
freq_unit = Unit of radio frequency (usually GHz) for plotting
time_offset_unit = Unit of time (usually s) for plotting
basename = Basename (string) for output plot (including path)
show_plot = Do you want to show the plot live? (True/False) (default = False)
'''
def plot_dedisp_subband_SNRvsDM(dedisp_ds,dedisp_timeseries,dedisp_times,freqs_array,t_cand,t_resol,trial_DMs,SNR,optimal_DM,offpulse_std_value,freq_unit='GHz',time_offset_unit='s',timeoffset_conversion_factor=1.0,basename='',show_plot=False):
    low_freq_limit = np.min(freqs_array)
    high_freq_limit = np.max(freqs_array)
    plot_name = basename+'_t'+'%.3f'% (t_cand)+'_dedispDS_SNRvsDM_freqs'+'%.2f'% (low_freq_limit)+'to'+'%.2f'% (high_freq_limit)+'.png'
    # Construct the gridspec framework for figure.
    fig = plt.figure(figsize=(6,8))
    #make outer gridspec
    outer = gridspec.GridSpec(2, 1, figure=fig,height_ratios = [1, 2])
    # Plot S/N vs. DM in the top gridspec.
    print('Plotting S/N vs. DM')
    gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = outer[0])
    ax1 = plt.subplot(gs1[0])
    ax1.plot(trial_DMs,SNR)
    ax1.set_xlabel('Trial DM (pc cm$^{-3}$)',fontsize=14)
    ax1.set_ylabel('S/N',fontsize=14)
    ax1.axvline(x=optimal_DM,label='DM$_{\mathrm{opt}}$',linestyle='--',color='r')
    ax1.legend(loc='best',prop={'size':12})
    # Construct nested gridspec within the bottom panel of outer gridspec.
    gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec = outer[1], height_ratios = [1, 3],hspace = 0)
    ax2 = plt.subplot(gs2[0])
    ax3 = plt.subplot(gs2[1],sharex=ax2)
    # Plot dedispersed time series in upper panel of bottom gridspec.
    print('Plotting dedispersed time series')
    # Convert dedispersed time series from integrated fluxes to S/N.
    dedisp_timeseries = dedisp_timeseries/offpulse_std_value
    # Convert times to offsets relative to pulse maximum.
    pulse_max_time = dedisp_times[np.argmax(dedisp_timeseries)]
    dedisp_time_offset = (dedisp_times -pulse_max_time)*timeoffset_conversion_factor # Dedispersed time offsets
    ax2.plot(dedisp_time_offset,dedisp_timeseries)
    ax2.set_ylabel('S/N',fontsize=14)
    # Plot the dedispersed dynamic spectrum in lower panel of bottom gridspec.
    print('Plotting dedispersed dynamic spectrum')
    ds_ext = [dedisp_time_offset[0],dedisp_time_offset[-1],freqs_array[0],freqs_array[-1]]
    ax3.imshow(dedisp_ds,origin='lower',interpolation='nearest',aspect='auto',extent=ds_ext)
    ax3.set_xlabel('Time offset ('+time_offset_unit+')',fontsize=14)
    ax3.set_ylabel('Radio frequency ('+freq_unit+')',fontsize=14)
    ax3.set_xlim(dedisp_time_offset[0],dedisp_time_offset[-1])
    plt.savefig(plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
# Combine plots of the entire dedispersed dynamic spectrum, dedispersed time series and S/N vs. DM plot.
'''
Inputs:
whole_ds = 2D array, DM = 0 pc/cc dynamic spectrum
times = 1D array of times (s) covered in non-dedispersed dynamic spectrum
dedisp_ds = 2D array, dedispersed dynamic spectrum
dedisp_timeseries = Sum of dedisp_ds over frequency
dedisp_times = 1D array of times(s) covered in dedisp_timeseries
freqs_array = 1D array of frequencies (GHz) covered in dedisp_ds
t_cand = Time (s) of occurrence of candidate
t_resol  = Time resolution (s) of data
trial_DMs = 1D array of trial DMs at which dedispersion must be performed
SNR = 1D array of signal-to-noise ratios at above trial DMs
optimal_DM = DM that maximized S/N
offpulse_std_value = Off-pulse standard deviation to normalize the dedispersed time series
flux_unit = Unit of flux density
freq_unit = Unit of radio frequency (usually GHz) for plotting
time_offset_unit = Unit of time (usually s) for plotting
timeoffset_conversion_factor = Conversion factor (float) from time units to time offset units. For example, s -> ms x 1.e3
basename = Basename (string) for output plot (including path)
show_plot = Do you want to show the plot live? (True/False) (default = False)
vmin = Min. color bar axis value for flux density (default = np.nanmin(whole_ds))
vmax = Max color bar axis value for flux density (default = np.nanmax(whole_ds))
log_colorbar = Do you want a log-spaced colorbar? (True/False) (default = False)
'''
def plot_dedisp_ds_SNRvsDM(whole_ds,times,dedisp_ds,dedisp_timeseries,dedisp_times,freqs_array,t_cand,t_resol,trial_DMs,SNR,optimal_DM,offpulse_std_value,flux_unit,freq_unit='GHz',time_offset_unit='s',timeoffset_conversion_factor=1.0,basename='',show_plot=False,vmin=None,vmax=None,log_colorbar=False):
    # Specify plot name.
    low_freq_limit = np.min(freqs_array)
    high_freq_limit = np.max(freqs_array)
    plot_name = basename+'_t'+'%.3f'% (t_cand)+'_dedispDS_SNRvsDM_freqs'+'%.2f'% (low_freq_limit)+'to'+'%.2f'% (high_freq_limit)+'.png'

    if (vmin is None):
        vmin = np.nanmin(whole_ds)
    elif isinstance(vmin, str):
        if ('median-' in vmin and 'sigma' in vmin):
            N = float(vmin.split('median-')[1].split('sigma')[0])
            vmin = np.nanmedian(whole_ds) - N*np.nanstd(whole_ds)
        elif ('mean-' in vmin and 'sigma' in vmin):
            N = float(vmin.split('mean-')[1].split('sigma')[0])
            vmin = np.nanmean(whole_ds) - N*np.nanstd(whole_ds)

    if (vmax is None):
        vmax = np.nanmax(whole_ds)
    elif isinstance(vmax, str):
        if ('median+' in vmax and 'sigma' in vmax):
            N = float(vmax.split('median+')[1].split('sigma')[0])
            vmax = np.nanmedian(whole_ds) + N*np.nanstd(whole_ds)
        elif ('mean+' in vmax and 'sigma' in vmax):
            N = float(vmax.split('mean+')[1].split('sigma')[0])
            vmax = np.nanmean(whole_ds) + N*np.nanstd(whole_ds)

    # Construct the gridspec framework for figure.
    fig = plt.figure(figsize=(8,14))
    #make outer gridspec
    outer = gridspec.GridSpec(3, 1, figure=fig,height_ratios = [1, 3, 2])

    # Plot S/N vs. DM in the top gridspec.
    print('Plotting S/N vs. DM')
    gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = outer[0])
    ax1 = plt.subplot(gs1[0])
    ax1.plot(trial_DMs,SNR)
    ax1.set_xlabel('Trial DM (pc cm$^{-3}$)',fontsize=16)
    ax1.set_ylabel('S/N',fontsize=16)
    ax1.axvline(x=optimal_DM,label='DM$_{\mathrm{opt}}$',linestyle='--',color='r')
    ax1.legend(loc='best',prop={'size':14})

    # Construct nested gridspec within the bottom panel of outer gridspec.
    gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec = outer[1], height_ratios = [1, 3],hspace = 0)
    ax2 = plt.subplot(gs2[0])
    ax3 = plt.subplot(gs2[1],sharex=ax2)

    # Determine time of pulse maximum.
    pulse_max_time = dedisp_times[np.argmax(dedisp_timeseries)]
    dedisp_time_offset = (dedisp_times -pulse_max_time)*timeoffset_conversion_factor
    times_offset = (times - pulse_max_time)*timeoffset_conversion_factor

    # Plot dedispersed time series in upper panel of bottom gridspec.
    print('Padding dedispersed data products with NaNs to the time span of the non-dedispersed dynamic spectrum')
    N_pad = whole_ds.shape[1] - np.size(dedisp_timeseries)
    dedisp_timeseries = np.append(dedisp_timeseries, np.ones(N_pad)*np.NaN)
    dedisp_ds = np.append(dedisp_ds, np.ones((len(whole_ds), N_pad))*np.NaN, axis=1)

    # Plot dedispersed time series in upper panel of bottom gridspec.
    print('Plotting dedispersed time series')
    dedisp_timeseries = dedisp_timeseries/offpulse_std_value # Convert dedispersed timeseries to S/N.
    ax2.plot(times_offset,dedisp_timeseries)
    ax2.set_ylabel('S/N',fontsize=16)

    # Plot the dedispersed dynamic spectrum in middle panel of bottom gridspec.
    print('Plotting dedispersed dynamic spectrum')
    ds_ext = [times_offset[0],times_offset[-1],freqs_array[0],freqs_array[-1]]
    if (log_colorbar==False):
        ax3.imshow(dedisp_ds,origin='lower',interpolation='nearest',aspect='auto',extent=ds_ext,vmin=vmin,vmax=vmax)
    else:
        ax3.imshow(dedisp_ds,origin='lower',interpolation='nearest',aspect='auto',extent=ds_ext,norm=LogNorm(vmin=vmin,vmax=vmax))
    ax3.set_ylabel('Radio frequency ('+freq_unit+')',fontsize=16)

    gs3 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = outer[2])
    ax4 = plt.subplot(gs3[0],sharex=ax2)
    # Plot non-dedispersed dynamic spectrum in lower panel of bottom gridspec.
    print('Plotting non-dedispersed dynamic spectrum')
    if (log_colorbar==False):
        im = ax4.imshow(whole_ds,origin='lower',interpolation='nearest',aspect='auto',extent=ds_ext,vmin=vmin,vmax=vmax)
    else:
        im = ax4.imshow(whole_ds,origin='lower',interpolation='nearest',aspect='auto',extent=ds_ext,norm=LogNorm(vmin=vmin,vmax=vmax))
    ax4.set_ylabel('Radio frequency ('+freq_unit+')',fontsize=16)
    ax4.set_xlabel('Time offset ('+time_offset_unit+') referenced to %.2f GHz'% (freqs_array[-1]),fontsize=16)
    ax4.set_xlim(times_offset[0],times_offset[-1])

    # Create space for colorbar on right side of plot.
    fig.subplots_adjust(right=0.80)
    ax4_left, ax4_bottom, ax4_width, ax4_height = ax4.get_position().get_points().flatten()
    ax2_left, ax2_bottom, ax2_width, ax2_height = ax2.get_position().get_points().flatten()
    cbar_ax = fig.add_axes([0.83, ax4_bottom, 0.03, ax2_bottom-ax4_bottom])
    h = fig.colorbar(im, cax=cbar_ax)
    h.set_label('Flux density (%s)'% (flux_unit),fontsize=16)

    plt.savefig(plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
# Plot a dedispersed dynamic spectrum together with the dedispersed time series.
'''
Inputs:
dedisp_ds = 2D array, dedispersed dynamic spectrum
dedisp_timeseries = 1D array, sum of dedisp_ds over frequency
dedisp_times = 1D array of times(s) covered in dedisp_timeseries
freqs_array = 1D array of frequencies (GHz) covered in dedisp_ds
t_cand = Time (s) of occurrence of candidate
t_resol  = Time resolution (s) of data
DM = Dispersion measure at which data were dedispersed
time_unit = Unit of time (usually s) for plotting
freq_unit = Unit of radio frequency (usually GHz) for plotting
flux_unit = Unit of flux density
basename = Basename (string) for output plot (including path)
show_plot = Do you want to show the plot live? (True/False) (default = False)
vmin = Min. color bar axis value for flux density (default = np.nanmin(ds))
vmax = Max color bar axis value for flux density (default = np.nanmax(ds))
log_colorbar = Do you want a log-spaced colorbar? (True/False) (default = False)
'''
def plot_dedispersed_ds(dedisp_ds,dedisp_timeseries,dedisp_times,freqs_array,t_cand,t_resol,DM,time_unit,freq_unit,flux_unit,basename,show_plot=False,vmin=None,vmax=None,log_colorbar=False):
    # Specify plot name.
    low_freq_limit = np.min(freqs_array)
    high_freq_limit = np.max(freqs_array)
    plot_name = basename+'_t'+'%.3f'% (t_cand)+'_dedispDS_DM_'+'%.1f'% (DM)+'_freqs'+'%.2f'% (low_freq_limit)+'to'+'%.2f'% (high_freq_limit)+'.png'
    if (vmin is None):
        vmin = np.nanmin(dedisp_ds)
    elif isinstance(vmin, str):
        if ('median-' in vmin and 'sigma' in vmin):
            N = float(vmin.split('median-')[1].split('sigma')[0])
            vmin = np.nanmedian(dedisp_ds) - N*np.nanstd(dedisp_ds)
        elif ('mean-' in vmin and 'sigma' in vmin):
            N = float(vmin.split('mean-')[1].split('sigma')[0])
            vmin = np.nanmean(dedisp_ds) - N*np.nanstd(dedisp_ds)

    if (vmax is None):
        vmax = np.nanmax(dedisp_ds)
    elif isinstance(vmax, str):
        if ('median+' in vmax and 'sigma' in vmax):
            N = float(vmax.split('median+')[1].split('sigma')[0])
            vmax = np.nanmedian(dedisp_ds) + N*np.nanstd(dedisp_ds)
        elif ('mean+' in vmax and 'sigma' in vmax):
            N = float(vmax.split('mean+')[1].split('sigma')[0])
            vmax = np.nanmean(dedisp_ds) + N*np.nanstd(dedisp_ds)

    print('Plotting dedispersed time series')
    ds_ext = [dedisp_times[0],dedisp_times[-1],freqs_array[0],freqs_array[-1]]
    fig,axes = plt.subplots(nrows=2,ncols=1,sharex=True,gridspec_kw={'height_ratios': [1, 2]})
    axes[0].plot(dedisp_times,dedisp_timeseries)
    axes[0].set_ylabel('Flux (%s)'% (flux_unit),fontsize=14)

    print('Plotting dedispersed dynamic spectrum')
    if (log_colorbar==False):
        im = axes[1].imshow(dedisp_ds,origin='lower',interpolation='nearest',aspect='auto',extent=ds_ext,vmin=vmin,vmax=vmax)
    else:
        im = axes[1].imshow(dedisp_ds,origin='lower',interpolation='nearest',aspect='auto',extent=ds_ext,norm=LogNorm(vmin=vmin,vmax=vmax))
    axes[1].set_xlabel('Time ('+time_unit+')',fontsize=14)
    axes[1].set_ylabel('Radio frequency ('+freq_unit+')',fontsize=14)
    axes[1].set_xlim((dedisp_times[0],dedisp_times[-1]))
    fig.subplots_adjust(hspace=0)

    # Create space for colorbar on right side of plot.
    fig.subplots_adjust(right=0.80)
    ax1_left, ax1_bottom, ax1_width, ax1_height = axes[1].get_position().get_points().flatten()
    ax0_left, ax0_bottom, ax0_width, ax0_height = axes[0].get_position().get_points().flatten()
    cbar_ax = fig.add_axes([0.83, ax1_bottom, 0.03, ax0_bottom-ax1_bottom])
    h = fig.colorbar(im, cax=cbar_ax)
    h.set_label('Flux density (%s)'% (flux_unit),fontsize=14)
    plt.savefig(plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
