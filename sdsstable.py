#!/usr/bin/env python
import numpy as np
import numpy.ma as ma
from astropy.io import ascii
from astropy.table import Table,join
from astropy import units as u
from astropy.coordinates import SkyCoord
#from astropy.visualization import quantity_support
#quantity_support()
#import matplotlib.pyplot as plt
#import pandas as pd

# Read in the SDSS standard star photemetry and astrometry catalogs 
# from Smith et al. 2002, AJ, 123, 2121
# and write them out as a combined table

photometrydata = 'tab08.dat.txt'
astrometrydata = 'datafile9.txt'

# create columns with the individual magnitudes and errors from the 
# columns with colors
t = Table.read(photometrydata,format='ipac')
t["u'"] = t["u'-g'"]+t["g'-r'"]+t["r'"]
t["rms u'"] = (t["rms r'"]**2 + t["rms g'-r'"]**2 + t["rms u'-g'"]**2)**.5
t["g'"] = t["r'"]+t["g'-r'"]
t["rms g'"] = (t["rms r'"]**2 + t["rms g'-r'"]**2)**.5
t["i'"] = t["r'"]-t["r'-i'"]
t["rms i'"] = (t["rms r'"]**2 + t["rms r'-i'"]**2)**.5
t["z'"] = t["i'"]-t["i'-z'"]
t["rms z'"] = (t["rms i'"]**2 + t["rms i'-z'"]**2)**.5

tastro = Table.read(astrometrydata,format='ipac')

# join the tables eliminating duplicate columns
tconcat= join(t,tastro,metadata_conflicts='warn')

#
# Make a coordinate pair from columns that look like coordinates :-)
#See http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord.guess_from_table
# an add a new column to table.
#tconcat['coordinates'] = SkyCoord.guess_from_table(t,unit=(u.hourangle,u.degree))

# Astropy IPAC format does not allow non-alphanumeric characters in column names, even with non-strict checking (DBMS=False),
# even though the actual IPAC standard allows them! http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html
formats = dict()
for c in tconcat.colnames:
    s = c.replace("'","").replace("-","_").replace(" ","_")
    tconcat[c].name = s
    # force a sensible output format for floats
    if tconcat[s].dtype == np.dtype('float64'): formats[s] = "%0.3f"
tconcat.write("ugriz_standards.tab",format='ipac',DBMS=False,formats=formats)


