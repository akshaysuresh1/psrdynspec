from .config import *
import matplotlib.gridspec as gridspec

# Plots saved to disk by default. Plots displayed on live window only if show_plot==True.
#########################################################################
# Plot single pulse candidates on DM vs. sigma plot.
'''
Inputs:
cand_DMs = DMs (pc/cc) of candidates
cand_sigma = Detection signficance of candidates
basename = Basename for plot
SAVE_DIR = Path to which plot must be saved
show_plot = Do you want to view the plot live? (True/False) (default = False)
'''
def scatter_DMvsSN(cand_DMs,cand_sigma,basename,SAVE_DIR,show_plot=False):
    plot_name = basename+'_scatterDMvsSN.png'
    plt.scatter(cand_DMs,cand_sigma,s=2)
    plt.xlabel('Dispersion measures (pc cm$^{-3}$)',fontsize=14)
    plt.ylabel('Detection S/N',fontsize=14)
    plt.savefig(SAVE_DIR+plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
#########################################################################
# Output a diagnostic DM vs. time plot for verifying uniquness of extracted single pulse candidates.
'''
Inputs:
cand_dedisp_times = Times (s) reported for all extracted candidates
cand_DMs = DMs (pc/cc) of all extracted candidates
unique_cand_dedisp_times = Times (s) reported for unique candidates
unique_cand_DMs = DMs (pc/cc) of unique candidates
cand_sigma = Detection significance of all extracted candidates
low_DM = Lowest DM (pc/cc) of interest
high_DM = Largest DM (pc/cc) of interest. Extracted single pulse candidates are identified at DMs between low_DM and high_DM (including endpoints).
time_margin = A candidate identified within this margin (s) of another candidate at slightly different DM is marked as a duplicate detection.
basename = Basename for output diagnostic plot
SAVE_DIR = Path to which diagnostic plot must be saved
show_plot = Do you want to show the diagnostic plot? (True/False) (default = False)
'''
def plot_diag_cands(cand_dedisp_times,cand_DMs,unique_cand_dedisp_times,unique_cand_DMs,cand_sigma,low_DM,high_DM,time_margin,basename,SAVE_DIR,show_plot=False):
    cand_dedisp_times = np.array(cand_dedisp_times,dtype=np.float64)
    cand_DMs = np.array(cand_DMs,dtype=np.float64)
    cand_sigma = np.array(cand_sigma,dtype=np.float64)

    # Plot the diagnostic DM vs. time plot.
    print('Plotting single pulse candidates on DM vs. time plane.')
    plot_name = basename+'_loDM'+str(low_DM)+'_hiDM'+str(high_DM)+'.png'
    plt.title('Minimum allowed time separation between unique candidates = %s ms'% (time_margin*1e3),fontsize=14)
    plt.plot(unique_cand_dedisp_times,unique_cand_DMs,marker='x',markersize=8,color='red',linestyle='None')
    plt.scatter(x=cand_dedisp_times,y=cand_DMs,s=(cand_sigma**2.)/5,marker='o',facecolor='None',edgecolor='k')
    plt.xlabel('Time (s)',fontsize=14)
    plt.ylabel('Dispersion Measure (pc cm$^{-3}$)',fontsize=14)
    plt.savefig(SAVE_DIR+plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
#########################################################################
# Reproduce plot of single pulse detections in the DM-time plane.
'''
Inputs:
cand_dedisp_times = 1D array of times (s) of single pulse candidates
cand_DMs = 1D array of DMs (pc/cc) of single pulse candidates
cand_sigma = 1D array of detection S/N of single pulse candidates
metadata = infodata object containing observation metadata
SAVE_DIR = Path to which output plot must be saved (d: working directory)
output_formats = List of file extensions for output plot (d: ['.png'])
show_plot = Do you want to show the plot live? (True/False) (default = False)
low_DM_cand = Min. DM (pc/cc) to plot (d: np.min(cand_DMs))
high_DM_cand = Max DM (pc/cc) to plot (d: np.max(cand_DMs))
'''
def plot_DMtime(cand_dedisp_times, cand_DMs, cand_sigma, metadata, SAVE_DIR='', output_formats=['.png'], show_plot=False, low_DM_cand=None, high_DM_cand=None):
    if low_DM_cand is None:
        low_DM_cand = np.min(cand_DMs)
    if high_DM_cand is None:
        high_DM_cand = np.min(cand_DMs)
    center_freq = metadata.lofreq + 0.5*(metadata.numchan-1)*metadata.chan_width # MHz

    plot_name = metadata.basename
    if '/' in plot_name:
        plot_name = plot_name.split('/')[-1]
    if 'DM' in plot_name:
        plot_name = plot_name.split('_DM')[0]

    print('Reconstructing single_pulse_search.py output')
    fig = plt.figure(figsize=(10,8))

    # Make outer gridspec,
    outer = gridspec.GridSpec(2, 1, figure=fig,height_ratios = [1, 2])
    # Plot DM-time plane in bottom gridspec.
    gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = outer[1])
    ax1 = plt.subplot(gs1[0])
    ax1.scatter(x=cand_dedisp_times,y=cand_DMs,s=(cand_sigma**2.)/5,marker='o',facecolor='None',edgecolor='k')
    ax1.set_ylim((low_DM_cand, high_DM_cand))
    ax1.set_xlabel('Time (s)',fontsize=14)
    ax1.set_ylabel('DM (pc cm$^{-3}$)',fontsize=14)
    # Split top gridspec into three subplots.
    gs0 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec = outer[0], width_ratios=[1,1,1],wspace=0.3)
    # Left subplot: Histogram of candidate S/N
    ax00 = plt.subplot(gs0[0])
    ax00.hist(cand_sigma, bins='auto', facecolor='None', edgecolor='k', histtype='step')
    ax00.set_xlim((np.min(cand_sigma),np.max(cand_sigma)))
    ax00.set_ylabel('No. of pulses per bin',fontsize=14)
    ax00.set_xlabel('S/N',fontsize=14)
    # Middle subplot: Histogram of candidate DMs
    ax01 = plt.subplot(gs0[1])
    ax01.hist(cand_DMs, bins='auto', facecolor='None', edgecolor='k', histtype='step')
    ax01.set_xlim((low_DM_cand, high_DM_cand))
    ax01.set_ylabel('No. of pulses per bin',fontsize=14)
    ax01.set_xlabel('DM (pc cm$^{-3}$)',fontsize=14)
    # Right subplot: Scatter plot candidate S/N vs. candidate DMs
    ax02 = plt.subplot(gs0[2])
    ax02.scatter(x=cand_DMs, y=cand_sigma, s=2, marker='o', facecolor='None', edgecolor='k')
    ax02.set_xlim((low_DM_cand, high_DM_cand))
    ax02.set_ylim((np.min(cand_sigma),np.max(cand_sigma)))
    ax02.set_xlabel('DM (pc cm$^{-3}$)',fontsize=14)
    ax02.set_ylabel('S/N', fontsize=14)
    fig.subplots_adjust(left=0.08, right=0.92, top=0.8, bottom=0.08)

    # Metadata
    # Column 0
    ax00.annotate('Source: %s'% (metadata.object),xy=(0.0, 1.5), xycoords='axes fraction',fontsize=14)
    ax00.annotate('Telescope: %s'% (metadata.telescope),xy=(0.0, 1.3), xycoords='axes fraction',fontsize=14)
    ax00.annotate('Instrument: %s'% (metadata.instrument),xy=(0.0, 1.1), xycoords='axes fraction',fontsize=14)
    # Column 1
    ax01.annotate('RA (J2000) = %s'% (metadata.RA),xy=(0.0, 1.5),xycoords='axes fraction',fontsize=14)
    if '-' in metadata.DEC:
        DEC_label = r'$-$' + metadata.DEC.split('-')[1]
    elif '+' in metadata.DEC:
        DEC_label = r'$+$' + metadata.DEC.split('+')[1]
    else:
        DEC_label = r'$+$' + metadata.DEC
    ax01.annotate('Dec (J2000) = %s'% (DEC_label),xy=(0.0, 1.3),xycoords='axes fraction',fontsize=14)
    if metadata.bary:
        MJD_label = r'MJD$_{\mathrm{bary}}$'
    else:
        MJD_label = r'MJD$_{\mathrm{topo}}$'
    ax01.annotate(MJD_label+' = %.12g'% (metadata.epoch),xy=(0.0, 1.1),xycoords='axes fraction',fontsize=14)
    # Column 2
    ax02.annotate('No. of samples = %d'% (metadata.N),xy=(0.0, 1.5), xycoords='axes fraction',fontsize=14)
    ax02.annotate('Sampling time = %.2f $\mu$s'% (metadata.dt*1e6),xy=(0.0, 1.3), xycoords='axes fraction',fontsize=14)
    ax02.annotate('Freq$_{\mathrm{ctr}}$ = %.1f MHz'% (center_freq),xy=(0.0, 1.1),xycoords='axes fraction',fontsize=14)

    # Save plot and display it live, if specified.
    for format in output_formats:
        plt.savefig(SAVE_DIR+plot_name+'_DMtime'+format)
    if show_plot:
        plt.show()
    else:
        plt.close()
#########################################################################
