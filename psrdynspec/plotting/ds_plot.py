from .config import *
from matplotlib.colors import LogNorm

# Plots saved to disk by default. Plots displayed on live window only if show_plot==True.
##########################################################################
# Plot a dynamic spectrum using imshow.
'''
ds = Dynamic spectrum to plot, 2D array, axes = [Frequency, Time]
times = 1D array of times (s) covered in dynamic spectrum
freqs_GHz = 1D array of radio frequencies (GHz) covered in dynamic spectrum
basename = Basename of output plot (including path)
show_plot = Do you want to show the plot live? (True/False) (default = False)
time_unit = String: Unit of time axis (s, ms, etc.)  (default = 's')
freq_unit = String: Unit of radio frequency (GHz, MHz, kHz, etc.) (default = 'GHz')
flux_unit = String: Unit of flux density measurement (default = 'arbitrary units')
vmin = Min. color bar axis value for flux density (default = np.nanmin(ds))
vmax = Max color bar axis value for flux density (default = np.nanmax(ds))
log_colorbar = Do you want a log-spaced colorbar? (True/False) (default = False)
cmap = Matplotlib color map (default = 'viridis')
mask_chans = List of channels to indicate as flagged channels (default = None)
'''
def plot_ds(ds,times,freqs_GHz,basename,show_plot=False,time_unit='s',freq_unit='GHz',flux_unit='arbitrary units',vmin=None,vmax=None,log_colorbar=False,cmap='viridis',mask_chans=None):
    plot_name = basename+'_t'+'%.2fto%.2f'% (times[0],times[-1])+'_freqs%.2fto%.2f'% (freqs_GHz[0],freqs_GHz[-1])+'.png'

    freq_resol = freqs_GHz[1] - freqs_GHz[0] # Frequency resolution

    if (vmin is None):
        vmin = np.nanmin(ds)
    elif isinstance(vmin, str):
        if ('median-' in vmin and 'sigma' in vmin):
            N = float(vmin.split('median-')[1].split('sigma')[0])
            vmin = np.nanmedian(ds) - N*np.nanstd(ds)
        elif ('mean-' in vmin and 'sigma' in vmin):
            N = float(vmin.split('mean-')[1].split('sigma')[0])
            vmin = np.nanmean(ds) - N*np.nanstd(ds)

    if (vmax is None):
        vmax = np.nanmax(ds)
    elif isinstance(vmax, str):
        if ('median+' in vmax and 'sigma' in vmax):
            N = float(vmax.split('median+')[1].split('sigma')[0])
            vmax = np.nanmedian(ds) + N*np.nanstd(ds)
        elif ('mean+' in vmax and 'sigma' in vmax):
            N = float(vmax.split('mean+')[1].split('sigma')[0])
            vmax = np.nanmean(ds) + N*np.nanstd(ds)

    # Set channels to mask in dynamic spectrum.
    if mask_chans is None:
        mask_chans = []

    print('Plotting dynamic spectrum...')
    ds_ext = [times[0], times[-1], freqs_GHz[0], freqs_GHz[-1]]
    if (log_colorbar==False):
        plt.imshow(ds,interpolation='nearest',origin='lower',aspect='auto',extent=ds_ext,cmap=cmap,vmin=vmin,vmax=vmax)
    else:
        plt.imshow(ds,interpolation='nearest',origin='lower',aspect='auto',extent=ds_ext,cmap=cmap,norm=LogNorm(vmin=vmin,vmax=vmax))
    # Indicate masked channels.
    mask_len = int(np.round(0.06*len(times)))
    for chan in mask_chans:
        plt.fill_between(times[:mask_len], freqs_GHz[chan]-0.5*freq_resol, freqs_GHz[chan]+0.5*freq_resol,color='salmon')
    h = plt.colorbar()
    h.set_label('Flux density ('+flux_unit+')',fontsize=14)
    plt.xlabel('Time ('+time_unit+')',fontsize=14)
    plt.ylabel('Radio frequency ('+freq_unit+')',fontsize=14)
    plt.savefig(plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
