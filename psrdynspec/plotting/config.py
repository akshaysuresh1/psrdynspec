# Set default plot settings.

import matplotlib as mpl
# Use Qt5Agg backend to enable interactive plotting (if needed).
mpl.use('Qt5Agg')
# Enable use of LaTeX labels in plots.
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Import required packages.
import matplotlib.pyplot as plt
plt.ioff() # Set interactive off by default.
import numpy as np
