# Set default plot settings.
from sys import platform
import matplotlib as mpl

# Set default backgrounds for plotting in screen/tmux without an X-window.
if platform=='darwin': # MacOS
    screen_backend = 'Agg'
    interactive_backend = 'MacOSX'
else:
    screen_backend = 'TkAgg'
    interactive_backend = 'Qt5Agg'

# If scripts are run from ipython, use interactive backend. Otherwise, use backends that do not require an X-window.
try:
    __IPYTHON__
except NameError:
    print("Not in IPython. Setting matplotlib backend to %s."% (screen_backend))
    mpl.use(screen_backend)
else:
    print("In IPython. Using %s backend for plotting"% (interactive_backend))
    mpl.use(interactive_backend)

# Enable use of LaTeX labels in plots.
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Import required packages.
import matplotlib.pyplot as plt
plt.ioff() # Set interactive off by default.
import numpy as np
