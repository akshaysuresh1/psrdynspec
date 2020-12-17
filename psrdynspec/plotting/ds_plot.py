from .config import *
from matplotlib.colors import LogNorm

# Plots saved to disk by default. Plots displayed on live window only if show_plot==True.
##########################################################################
# Plot a dynamic spectrum using imshow.
'''
ds = Dynamic spectrum to plot, 2D array, axes = [Frequency, Time]
t_start = Start time of DS
t_stop = Stop time of DS
freq_low = Frequency at lower edge of DS bandwidth
freq_high = Frequency at upper edge of DS bandwidth
time_unit = String: Unit of time axis (s, ms, etc.)
freq_unit = String: Unit of radio frequency (GHz, MHz, kHz, etc.)
flux_unit = String: Unit of flux density measurement
basename = Basename of output plot (including path)
show_plot = Do you want to show the plot live? (True/False) (default = False)
vmin = Min. color bar axis value for flux density (default = np.nanmin(ds))
vmax = Max color bar axis value for flux density (default = np.nanmax(ds))
log_colorbar = Do you want a log-spaced colorbar? (True/False) (default = False)
cmap = Matplotlib color map (default = 'viridis')
mask_chans = List of channels to indicate as flagged channels (default = None)
'''
def plot_ds(ds,t_start,t_stop,freq_low,freq_high,time_unit,freq_unit,flux_unit,basename,show_plot=False,vmin=None,vmax=None,log_colorbar=False,cmap='viridis',mask_chans=None):
    plot_name = basename+'_t'+'%.3fto%.3f'% (t_start,t_stop)+'_freqs%.2fto%.2f'% (freq_low,freq_high)+'.png'
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
    if (log_colorbar==False):
        plt.imshow(ds,interpolation='nearest',origin='lower',aspect='auto',extent=[t_start,t_stop,freq_low,freq_high],cmap=cmap,vmin=vmin,vmax=vmax)
    else:
        plt.imshow(ds,interpolation='nearest',origin='lower',aspect='auto',extent=[t_start,t_stop,freq_low,freq_high],cmap=cmap,norm=LogNorm(vmin=vmin,vmax=vmax))
    for chan in mask_chans:
        plt.axhline(y=freqs_GHz[chan],xmin=0., xmax=0.03, linestyle='-', color='salmon')
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
