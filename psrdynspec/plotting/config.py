# Set default plot settings.
from sys import platform
import matplotlib as mpl

if platform=='darwin': # MacOS
    noninteractive_backend = 'Agg' # Set default backgrounds for non-interactive plotting on screen/tmux.
    ipython_shell_backend = 'MacOSX'
else:
    noninteractive_backend = 'Agg' # Set default backgrounds for non-interactive plotting on screen/tmux.
    ipython_shell_backend = 'Qt5Agg'

# If scripts are run from ipython, use interactive backend. Otherwise, use non-interactive backends.
# Use 'nbAgg' backend when run within a Jupyter notebook.
try:
    __IPYTHON__
except NameError:
    print("Enabling non-interactive plotting. Setting matplotlib backend to %s."% (screen_backend))
    mpl.use(noninteractive_backend)
else:
    ipython_class = get_ipython().__class__.__name__
    if ipython_class=='TerminalInteractiveShell':
        print("In Ipython shell. Using %s backend for plotting"% (ipython_shell_backend))
        mpl.use(ipython_shell_backend)
    elif ipython_class=='ZMQInteractiveShell':
        print("In Jupyter notebook. Using ngAgg backend for plotting")
        mpl.use('nbAgg')
    else:
        print("Enabling non-interactive plotting. Setting matplotlib backend to %s."% (screen_backend))
        mpl.use(noninteractive_backend)

# Enable use of LaTeX labels in plots.
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Import required packages.
import matplotlib.pyplot as plt
plt.ioff() # Set interactive off by default.
import numpy as np
