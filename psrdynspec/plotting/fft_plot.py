from .config import *
import matplotlib.gridspec as gridspec

# Plots saved to disk by default. Plots displayed on live window only if show_plot==True.
##########################################################################
# Plot FFT of 1D time series.
'''
Inputs:
fourierfreqs = 1D array of Fourier frequencies
powerspec = 1D array of power spectra values
fourierfreq_units = Units of Fourier frequency (Hz, kHz, MHz, etc..)
powerspec_units = Units of the power spectrum
basename = Basename of output plot
SAVE_DIR = Path to which plot must be saved
show_plot = Do you want to show the plot? (True/False) (default = False)
axes_format = 'linear' / 'loglog' / 'semilogx'/'semilogy' (default = 'linear')
positive_freqs_only = Plot only the positive frequencies? (True/False) (default = False)
'''
def plot_FFT1d(fourierfreqs,powerspec,fourierfreq_units,powerspec_units,basename,SAVE_DIR,show_plot=False,axes_format='linear',positive_freqs_only=False):
    if (positive_freqs_only==True):
        pos_indices = np.where(fourierfreqs>0)[0]
        fourierfreqs = fourierfreqs[pos_indices]
        powerspec = powerspec[pos_indices]

    print('Plotting power spectrum...')
    plot_name = basename+'_FFT.png'
    plt.figure()
    if (axes_format=='semilogx'):
        plt.semilogx(fourierfreqs,powerspec)
    elif (axes_format=='semilogy'):
        plt.semilogy(fourierfreqs,powerspec)
    elif (axes_format=='loglog'):
        plt.loglog(fourierfreqs,powerspec)
    else:
        plt.plot(fourierfreqs,powerspec)
    plt.xlabel('Frequency (%s)'% (fourierfreq_units),fontsize=14)
    plt.ylabel('Spectral power density (%s)'% (powerspec_units),fontsize=14)
    plt.savefig(SAVE_DIR+plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
# Plot 1D time series and its FFT over different ranges of fourier frequency.
'''
Inputs:
times = 1D array of times (s)
timeseries = 1D array of flux value
fourier_freqs = 1D array of fourier frequencies (Hz) returned by FFT
power_spectrum = 1D array of power spectrum values returned by FFT
max_fourierfreq_plot = Max. fourier frequency (Hz) to be shown in plot
timeseries_unit = Unit of timeseries flux
powerspec_unit = Unit of power spectrum measurement
DM = Dispersion measure (pc/cc) used to dedisperse the time series
radiofreq_annotation = String to indicate radio frequencies integrated to produce plotted timeseries
special_fourierfreq = 1D array of specific frequencies (Hz) to indicate via vertical dashed lines in FFT plot
basename = Basename of output plot
SAVE_DIR = Path to which plot must be saved
show_plot = Do you want to show the plot? (True/False) (default = False)
plot_format = Format of output plot E.g., '.png', '.eps', etc.
'''
def fft_gridplot(times, timeseries, fourier_freqs, power_spectrum, max_fourierfreq_plot, timeseries_unit, powerspec_unit, DM, radiofreq_annotation, special_fourierfreq, basename, SAVE_DIR, show_plot, plot_format):
    plot_name = basename+'_FFTgrid'+plot_format
    # Construct the gridspec framework for figure.
    fig = plt.figure(figsize=(8,8))
    #make outer gridspec
    outer = gridspec.GridSpec(2, 1, figure=fig,height_ratios = [1, 4])

    # Plot timeseries in top gridspec.
    gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = outer[0])
    ax1 = plt.subplot(gs1[0])
    ax1.plot(times, timeseries,color='k')
    if np.ma.is_masked(timeseries):
        plot_name = plot_name = basename+'_noisemasked_FFTgrid'+plot_format
        ax1.axhline(y=np.min(timeseries),linestyle='--',color='orange')    
    ax1.annotate(radiofreq_annotation, xy=(0.05,0.8), xycoords='axes fraction',fontsize=14)
    ax1.annotate('DM = %.1f pc cm$^{-3}$'% (DM), xy=(0.7,0.8), xycoords='axes fraction',fontsize=14)
    ax1.set_xlabel('Time (s)', fontsize=14)
    ax1.set_ylabel('Flux (%s)'% (timeseries_unit), fontsize=14)
    ax1.set_ylim((np.min(timeseries),2*np.max(timeseries)))

    # Split bottom gridspec into 4 sub-panels for plotting FFT over different ranges of fourier frequency.
    gs2 = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec = outer[1],hspace=0.3)
    # Full frequency range in top sub-panel.
    ax20 = plt.subplot(gs2[0])
    ax20.plot(fourier_freqs, power_spectrum,color='k')
    ax20.set_xlim((0,max_fourierfreq_plot))
    # Show FFT over [0, 0.5*max_fourierfreq_plot] in second sub-panel.
    ax21 = plt.subplot(gs2[1])
    ax21.plot(fourier_freqs, power_spectrum,color='k')
    ax21.set_xlim((0,0.5*max_fourierfreq_plot))
    # Show FFT over [0, 0.25*max_fourierfreq_plot] in second sub-panel.
    ax22 = plt.subplot(gs2[2])
    ax22.plot(fourier_freqs, power_spectrum,color='k')
    ax22.set_xlim((0,0.25*max_fourierfreq_plot))
    # Show FFT over 0-10 Hz in second sub-panel.
    ax23 = plt.subplot(gs2[3])
    ax23.plot(fourier_freqs, power_spectrum,color='k')
    ax23.set_xlim((0,10))
    ax23.set_xlabel('Fourier frequency (Hz)',fontsize=14)

    # Indicate special fourier frequencies in all FFT panels.
    for i in range(len(special_fourierfreq)):
        freq = special_fourierfreq[i]
        ax20.axvline(x=freq,linestyle='--',color='grey')
        ax21.axvline(x=freq,linestyle='--',color='grey')
        ax22.axvline(x=freq,linestyle='--',color='grey')
        ax23.axvline(x=freq,linestyle='--',color='grey')

    fig.text(0.03, 0.4, 'Power spectrum (%s)'% (powerspec_unit), va='center', rotation='vertical',fontsize=14)
    fig.subplots_adjust(left=0.1, right=0.92, bottom=0.1,top=0.9)
    plt.savefig(SAVE_DIR+plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
