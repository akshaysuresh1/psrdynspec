# Define 'infodata' class to read and write information from PRESTO .inf files.

'''
The following script is adapted from PRESTO (https://github.com/scottransom/presto/blob/master/python/presto/infodata.py) under the following license.

GNU General Public License
Version 2, June 1991

Copyright (C) 1998-2016 Scott M. Ransom <sransom@nrao.edu>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
'''
##############################################################################
from builtins import object

class infodata(object):
    # Read a PRESTO .inf file and initialize an infodata object to organize its contents.
    def __init__(self, filenm):
        self.breaks = 0
        for line in open(filenm):
            if line.startswith(" Data file name"):
                self.basename = line.split("=")[-1].strip()
                continue
            if line.startswith(" Telescope"):
                self.telescope = line.split("=")[-1].strip()
                continue
            if line.startswith(" Instrument"):
                self.instrument = line.split("=")[-1].strip()
                continue
            if line.startswith(" Object being observed"):
                self.object = line.split("=")[-1].strip()
                continue
            if line.startswith(" J2000 Right Ascension"):
                self.RA = line.split("=")[-1].strip()
                continue
            if line.startswith(" J2000 Declination"):
                self.DEC = line.split("=")[-1].strip()
                continue
            if line.startswith(" Data observed by"):
                self.observer = line.split("=")[-1].strip()
                continue
            if line.startswith(" Epoch"):
                self.epoch = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Barycentered?"):
                self.bary = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Number of bins"):
                self.N = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Width of each time series bin"):
                self.dt = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Any breaks in the data?"):
                self.breaks = int(line.split("=")[-1].strip())
                if self.breaks:
                    self.onoff = []
                continue
            if line.startswith(" On/Off bin pair"):
                vals = line.split("=")[-1].strip().split(",")
                self.onoff.append((int(vals[0]), int(vals[1])))
                continue
            if line.startswith(" Type of observation"):
                self.waveband = line.split("=")[-1].strip()
                continue
            if line.startswith(" Beam diameter"):
                self.beam_diam = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Dispersion measure"):
                self.DM = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Central freq of low channel"):
                self.lofreq = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Total bandwidth"):
                self.BW = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Number of channels"):
                self.numchan = int(line.split("=")[-1].strip())
                continue
            if line.startswith(" Channel bandwidth"):
                self.chan_width = float(line.split("=")[-1].strip())
                continue
            if line.startswith(" Data analyzed by"):
                self.analyzer = line.split("=")[-1].strip()
                continue

    # Write infodata object to a .inf file readable by PRESTO.
    def to_file(self, inffn, notes=None):
        if not inffn.endswith(".inf"):
            raise ValueError("PRESTO info files must end with '.inf'. "
                             "Got: %s" % inffn)
        with open(inffn, 'w') as ff:
            if hasattr(self, 'basename'):
                ff.write(" Data file name without suffix          =  %s\n" %
                         self.basename)
            if hasattr(self, 'telescope'):
                ff.write(" Telescope used                         =  %s\n" %
                         self.telescope)
            if hasattr(self, 'instrument'):
                ff.write(" Instrument used                        =  %s\n" %
                         self.instrument)
            if hasattr(self, 'object'):
                ff.write(" Object being observed                  =  %s\n" %
                         self.object)
            if hasattr(self, 'RA'):
                ff.write(" J2000 Right Ascension (hh:mm:ss.ssss)  =  %s\n" %
                         self.RA)
            if hasattr(self, 'DEC'):
                ff.write(" J2000 Declination     (dd:mm:ss.ssss)  =  %s\n" %
                         self.DEC)
            if hasattr(self, 'observer'):
                ff.write(" Data observed by                       =  %s\n" %
                         self.observer)
            if hasattr(self, 'epoch'):
                ff.write(" Epoch of observation (MJD)             =  %05.15f\n" %
                         self.epoch)
            if hasattr(self, 'bary'):
                ff.write(" Barycentered?           (1=yes, 0=no)  =  %d\n" %
                         self.bary)
            if hasattr(self, 'N'):
                ff.write(" Number of bins in the time series      =  %-11.0f\n" %
                         self.N)
            if hasattr(self, 'dt'):
                ff.write(" Width of each time series bin (sec)    =  %.15g\n" %
                         self.dt)
            if hasattr(self, 'breaks') and self.breaks:
                ff.write(" Any breaks in the data? (1 yes, 0 no)  =  1\n")
                if hasattr(self, 'onoff'):
                    for ii, (on, off) in enumerate(self.onoff, 1):
                        ff.write(" On/Off bin pair #%3d                   =  %-11.0f, %-11.0f\n" %
                                 (ii, on, off))
            else:
                ff.write(" Any breaks in the data? (1 yes, 0 no)  =  0\n")
            if hasattr(self, 'DM'):
                ff.write(" Dispersion measure (cm-3 pc)           =  %.12g\n" %
                         self.DM)
            if hasattr(self, 'lofreq'):
                ff.write(" Central freq of low channel (Mhz)      =  %.12g\n" %
                         self.lofreq)
            if hasattr(self, 'BW'):
                ff.write(" Total bandwidth (Mhz)                  =  %.12g\n" %
                         self.BW)
            if hasattr(self, 'numchan'):
                ff.write(" Number of channels                     =  %d\n" %
                         self.numchan)
            if hasattr(self, 'chan_width'):
                ff.write(" Channel bandwidth (Mhz)                =  %.12g\n" %
                         self.chan_width)
            if hasattr(self, 'analyzer'):
                ff.write(" Data analyzed by                       =  %s\n" %
                         self.analyzer)
            if hasattr(self, 'deorbited'):
                ff.write(" Orbit removed?          (1=yes, 0=no)  =  %d\n" %
                         self.deorbited)
            ff.write(" Any additional notes:\n")
            if notes is not None:
                ff.write("    %s\n" % notes.strip())

    # Print attributes of an infodata object.
    def __str__(self):
        # Check if the timeseries contains any breaks.
        if hasattr(self, 'onoff'):
            breaks_line = ["Any breaks in the data? (1 yes, 0 no)  =  1"]
            for ii, (on, off) in enumerate(self.onoff, 1):
                onoff_description = "On/Off bin pair #%3d                   =  %-11.0f, %-11.0f" % (ii, on, off)
                breaks_line.append(onoff_description)
            breaks_line = '\n'.join(breaks_line)
        else:
            breaks_line = "Any breaks in the data? (1 yes, 0 no)  =  0"
        # Print output lines following the same structure as the .inf file content.
        lines = [
        "Data file name without suffix          =  %s"% (self.basename),
        "Telescope used                         =  %s"% (self.telescope),
        "Instrument used                        =  %s"% (self.instrument),
        "Object being observed                  =  %s"% (self.object),
        "J2000 Right Ascension (hh:mm:ss.ssss)  =  %s"% (self.RA),
        "J2000 Declination     (dd:mm:ss.ssss)  =  %s"% (self.DEC),
        "Data observed by                       =  %s"% (self.observer),
        "Epoch of observation (MJD)             =  %05.15f"% (self.epoch),
        "Barycentered?           (1=yes, 0=no)  =  %d"% (self.bary),
        "Number of bins in the time series      =  %-11.0f"% (self.N),
        "Width of each time series bin (usec)   =  %.15g"% (self.dt*1e6),
        breaks_line,
        "Type of observation (EM band)          =  %s"% (self.waveband),
        "Beam diameter (arcsec)                 =  %s"% (self.beam_diam),
        "Dispersion measure (pc/cc)             =  %.12g"% (self.DM),
        "Central freq of low channel (MHz)      =  %.12g"% (self.lofreq),
        "Total bandwidth (MHz)                  =  %.12g"% (self.BW),
        "Number of channels                     =  %d"% (self.numchan),
        "Channel bandwidth (MHz)                =  %.12g"% (self.chan_width),
        "Data analyzed by                       =  %s"% (self.analyzer)
        ]
        output = '\n'.join(lines)
        return output

    # Object representation
    def __repr__(self):
        return str(self)

##############################################################################
