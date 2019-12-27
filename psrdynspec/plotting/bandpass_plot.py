from .config import *

# Plots saved to disk by default. Plots displayed on live window only if show_plot==True.
##########################################################################
# Plot bandpass shape as a function of radio frequency.
'''
Inputs:
freqs_GHz = 1D array of radio frequencies (usually GHz)
bandpass = 1D array of bandpass fluxes
freq_unit = String: Unit of radio frequency (GHz, MHz, kHz, etc.)
flux_unit = String: Unit of flux density measurement
basename = Basename of output plot (string)
SAVE_DIR = Path to which output plots must be saved (string)
show_plot = Do you want to show the plot live? (True/False) (default = False)
'''
def plot_bandpass(freqs_GHz,bandpass,freq_unit,flux_unit,basename,SAVE_DIR,show_plot=False):
    plot_name = basename+'_bandpass.png'
    print('Plotting bandpass shape as a function of radio frequency')
    plt.plot(freqs_GHz,bandpass)
    plt.xlabel('Radio frequency ('+freq_unit+')',fontsize=14)
    plt.ylabel('Flux density ('+flux_unit+')',fontsize=14)
    plt.savefig(SAVE_DIR+plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
# Plot the bandpass shape for selected subbands.
'''
Inputs:
freqs_GHz = 1D array of radio frequencies (usually GHz)
bandpass = 1D array of bandpass fluxes over entire frequency range
freqs_subband = List of 1D arrays of frequency values spanned by subbands, len(freqs_subband) = No. of subbands
bandpass_fit_subband = List of 1D arrays of bandpass fits performed independently over each subband
freq_unit = String: Unit of radio frequency (GHz, MHz, kHz, etc.)
flux_unit = String: Unit of flux density measurement
basename = Basename of output plot (string)
SAVE_DIR = Path to which output plots must be saved (string)
show_plot = Do you want to show the plot live? (True/False) (default = False)
'''
def plot_bandpass_subbands(freqs_GHz,median_bp,freqs_subband,bandpass_fit_subband,freq_unit,flux_unit,basename,SAVE_DIR,show_plot=False):
    plot_name = basename+'_bandpass.png'
    print('Plotting median bandpass shape as a function of radio frequency')
    plt.plot(freqs_GHz,median_bp,color='b',alpha=0.4)
    for i in range(len(freqs_subband)):
        print('Plotting bandpass fit for subband %d: %.2f - %.2f %s'% (i+1,np.min(freqs_subband[i]),np.max(freqs_subband[i]),freq_unit))
        plt.axvspan(xmin=freqs_subband[i][0],xmax=freqs_subband[i][-1],color='grey',alpha=0.1)
        plt.plot(freqs_subband[i],bandpass_fit_subband[i],color='r')
    plt.xlabel('Radio frequency ('+freq_unit+')',fontsize=14)
    plt.ylabel('Flux density ('+flux_unit+')',fontsize=14)
    plt.savefig(SAVE_DIR+plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
