from __future__ import absolute_import
from pkg_resources import get_distribution, DistributionNotFound

try:
    from .io.parse_sp import read_singlepulse, gen_singlepulse
    from .io.read_config import read_config
    from .io.read_fil import load_fil_data
    from .io.read_header import Header
    from .io.read_psrfits import load_psrfits_data

except ImportError:
    print("Warning: Cannot import main I/O utilities.")

try:
    from .modules import dedisperse, ds_systematics, fft, filters1d, filters2d, fold, subbands
except ImportError:
    print("Warning: Cannot import main processing utilities.")

try:
    __version__ = get_distribution('psrdynspec').version
except DistributionNotFound:
    __version__ = '0.0.0 - please install via setup.py'
