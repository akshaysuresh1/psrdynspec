from .config import *

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
