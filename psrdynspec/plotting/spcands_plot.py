from .config import *

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
