from .config import *

# Plots saved to disk by default. Plots displayed on live window only if show_plot==True.
##########################################################################
# Plot 1D data, calculated trend and detrended data.
'''
Inputs:
axis_values = 1D array of axis values
data_1d = 1D array of data values
trend =  1D array of values representing calculated trend
detrended_data = 1D array of detrended data values
basename = Basename of output plot
xlabel = x-label of plot
ylabel = y-label of plot
SAVE_DIR = Path to which plot must be saved
show_plot = Do you want to show the plot live? (True/False) (default = False)
'''
def plot_detrended_data1D(axis_values,data_1d,trend,detrended_data,basename,xlabel,ylabel,SAVE_DIR,show_plot=False):
    plot_name = basename+'_detrend.png'
    print('Plotting data and calculated trend...')
    fig,axes = plt.subplots(nrows=2,ncols=1,sharex=True,gridspec_kw={'height_ratios': [1, 1]},figsize=((6,8)))
    axes[0].plot(axis_values,data_1d,label='Data')
    axes[0].plot(axis_values,trend,label='Subtracted trend')
    axes[0].set_ylabel(ylabel,fontsize=14)
    axes[0].legend(loc='best',prop={'size':10})

    print('Plotting detrended data...')
    axes[1].plot(axis_values,detrended_data,label='Detrended data')
    axes[1].set_xlabel(xlabel,fontsize=14)
    axes[1].set_ylabel(ylabel,fontsize=14)
    axes[1].legend(loc='best',prop={'size':10})
    axes[1].set_xlim((np.floor(np.min(axis_values)),np.ceil(np.max(axis_values))))
    fig.subplots_adjust(hspace=0)
    plt.savefig(SAVE_DIR+plot_name)
    if (show_plot==True):
        plt.show()
    else:
        plt.close()
##########################################################################
