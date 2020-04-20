from .config import *

# Plots saved to disk by default. Plots displayed on live window only if show_plot==True.
##########################################################################
# Assign plotting labels to different metrics.
'''
Inputs:
metric = Metric to maximize for determining optimal folding period ('reducedchisquare', 'profmax', 'profSNR')
'''
def assign_metlabel(metric):
    if metric=='profmax':
        label = 'Profile maximum'
    elif metric=='reducedchisquare':
        label = r'$\chi_{\mathrm{r}}^2$'
    elif metric=='profSNR':
        label = 'Profile S/N'
    return label
###########################################################################
# Plot metric values vs. trial folding periods.
'''
Inputs:
trial_periods = 1D array of trial folding periods (s)
metric_values = 1D array of metric values at above folding periods
metric = Metric to maximize for determining optimal folding period ('reducedchisquare', 'profmax', 'profSNR')
basename = Basename of output plot
SAVE_DIR = Path to which plot must be saved
show_plot = Do you want to show the plot live? (True/False) (default = False)
'''
def plot_metricvsperiod(trial_periods,metric_values,metric,basename,SAVE_DIR,show_plot=False):
    plot_name = basename+'_%s_vs_periods.png'% (metric)
    met_label = assign_metlabel(metric)
    print('Plotting metric values vs. trial folding periods...')
    plt.figure()
    plt.plot(trial_periods,metric_values)
    plt.xlabel('Trial folding period (s)',fontsize=14)
    plt.ylabel(met_label,fontsize=14)
    plt.savefig(SAVE_DIR+plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
# Plot a folded pulse profile.
'''
Inputs:
phasebins = 1D array of phase values
profile = Average pulse profile (1D array)
pfold = Folding period (s)
basename = Basename of output plot
SAVE_DIR = Path to which plot must be saved
show_plot = Do you want to show the plot live? (True/False)  (default = False)
'''
def plot_avgpulseprofile(phasebins,profile,pfold,basename,SAVE_DIR,show_plot=False):
    print('Plotting average pulse profile for P = %.5f s'% (pfold))
    plot_name = basename+'_avgprofile_P%.5f'% (pfold)+'.png'
    plt.figure()
    plt.plot(phasebins,profile)
    plt.title('Average pulse profile obtained using P = %.5f s'% (pfold),fontsize=14)
    plt.xlabel('Phase',fontsize=14)
    plt.ylabel('Flux (arbitrary units)',fontsize=14)
    plt.savefig(SAVE_DIR+plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
# Combine subplots of metric values vs. trial periods and folded pulse profile at optimal period in a single figure.
'''
Inputs:
trial_periods = 1D array of trial folding periods (s)
metric_values = 1D array of metric values at above folding periods
metric = Metric to maximize for determining optimal folding period ('reducedchisquare', 'profmax', 'profSNR')
phasebins = 1D array of phase values
profile = Average pulse profile (1D array) obtained by folding time series at the best period
best_period = Folding period (s) that maximizes the chosen metric.
basename = Basename of output plot
SAVE_DIR = Path to which plot must be saved
show_plot = Do you want to show the plot live? (True/False)  (default = False)
'''
def subplots_metric_profile(trial_periods,metric_values,metric,phasebins,profile,best_period,basename,SAVE_DIR,show_plot=False):
    plot_name = basename+'_%s_bestprofile.png'% (metric)
    met_label = assign_metlabel(metric)
    fig,axes = plt.subplots(nrows=2,ncols=1,gridspec_kw={'height_ratios': [1, 1.5]},figsize=(9,6.5))
    axes[0].plot(trial_periods,metric_values, marker='o', markersize=2, alpha=0.5,color='k')
    #axes[0].axvline(x=best_period,label='P$_{\mathrm{opt}}$ = %.5f s'% (best_period),linestyle='--',color='darkorange')
    axes[0].set_xlim(np.min(periods), np.max(periods))
    axes[0].set_xlabel('Trial folding period (s)',fontsize=14)
    axes[0].set_ylabel(met_label,fontsize=14)
    #axes[0].legend(loc='best',prop={'size':14})
    print('Plotting folded pulse profile at best period.')
    axes[1].plot(phasebins,profile,'-k')
    axes[1].annotate('P$_{\mathrm{opt}}$ = %.5f s'% (best_period), xycoords='axes fraction',xy=(0.8,0.91),fontsize=14)
    axes[1].set_xlabel('Phase',fontsize=14)
    axes[1].set_ylabel('Flux (arbitrary units)',fontsize=14)
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9,bottom=0.1,hspace=0.3)
    plt.savefig(SAVE_DIR+plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
# Plot folded pulse profile for each rotation as a function of phase.
'''
Inputs:
profile_rotations = Folded pulse profile for every rotation (2D array = (rotation number, phase bin))
counts_perrot_phibin = Counts per rotation per phase bin (2D array = (rotation number, phase bin))
phibins = Phase bins (1D array)
basename = Basename of output plot
SAVE_DIR = Path to which plot must be saved
show_plot = Do you want to show the plot live? (True/False) (default = False)
low_phase_limit = Lower limit to phase axis for plotting (default = 0.0)
high_phase_limit = Upper limit to phase axis for plotting (default = 1.0)
rot_spacing = Spacing between consecutive rotations on y-axis of plot (default = 1.0)
normalization = Pulse normalization factor (float) to aid plotting (default = 'quarterrotmax', i.e., (0.25*(Total no.of rotations)/ global pulse max) )
'''
def plot_foldedprofile_rotations(profile_rotations,counts_perrot_phibin,phibins,basename,SAVE_DIR,show_plot = False,low_phase_limit=0.0,high_phase_limit=1.0,rot_spacing = 1.0, normalization = 'quarterrotmax'):
    N_rotations = len(profile_rotations)
    plot_name = basename+'_folded_rotations.png'

    # Scale pulse profiles across all rotations by this factor to aid plotting.
    if (normalization=='quarterrotmax'):
        rot_scaling_factor = 0.25*N_rotations/np.nanmax(profile_rotations)
    else:
        rot_scaling_factor = normalization

    fig,axes = plt.subplots(nrows=2,ncols=1,sharex=True,gridspec_kw={'height_ratios': [1, 2]},figsize=(6,7.2))
    integrated_profile = np.nansum(profile_rotations*counts_perrot_phibin,axis=0)/np.nansum(counts_perrot_phibin,axis=0)
    axes[0].plot(phibins,integrated_profile,'-k')
    axes[0].set_ylabel('Flux (arbitrary units)',fontsize=14)
    profile_rotations = profile_rotations*rot_scaling_factor
    for i in range(N_rotations):
        axes[1].plot(phibins,rot_spacing*i+profile_rotations[i],'-k',linewidth=1)
    axes[1].set_xlabel('Phase',fontsize=14)
    axes[1].set_ylabel('Rotations',fontsize=14)
    axes[1].set_xlim((low_phase_limit,high_phase_limit))
    axes[1].set_ylim((0,np.ceil(N_rotations*rot_spacing/10)*10))
    # Adjust y-axis ticklabels if more spacing between profiles per rotation is desired.
    if (rot_spacing!=1):
        yticklabels = np.arange(0,np.ceil(N_rotations/10)*10+1,10)
        yticks = yticklabels*rot_spacing
        axes[1].set_yticks(yticks)
        axes[1].set_yticklabels(yticklabels)
    fig.subplots_adjust(hspace=0.075)
    plt.savefig(SAVE_DIR+plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
