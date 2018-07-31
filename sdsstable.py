#!/usr/bin/env python
import numpy as np
import numpy.ma as ma
from astropy.io import ascii
from astropy.table import Table,join,Column
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import re
import filtermanage as fm

#from astropy.visualization import quantity_support
#quantity_support()
#import matplotlib.pyplot as plt
#import pandas as pd

# Read in the SDSS standard star photemetry and astrometry catalogs 
# from Smith et al. 2002, AJ, 123, 2121
# and write them out as a combined table

photometrydata = 'tab08.dat.txt'
astrometrydata = 'datafile9.txt'
# this table has UVBRIJHK and distances
moredata       = 'sdss_standards_UBVRIJHK_distance.xml'

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

# can't join moretable because the sources have alternate names.
# So create a column called StarName in moretable to allow a join.
# Also make that an actual string not bytes. I hate bytes
moretable = Table.read(moredata,format='votable')
s = []
for i in range(len(moretable)):
    c = moretable['MAIN_ID'][i]
    s.append(re.sub('\s+',' ',c.decode("utf-8")))
moretable.add_column(Column(data = s,name='StarName'))

tconcat.add_index('StarName')
#moretable.add_index('StarName')

#--------------
# playing with pandas to match entries, not needed, as coords below
# is more robust
#df1 = tconcat.to_pandas()
#df2 = moretable.to_pandas()
#
#match = 0
#nomatch = 0
#for name in moretable['StarName']:
#   if True:
#    if len(df1.loc[df1['StarName'] == name]) != 0:
#        print("Found match: %s"%name)
#        match += 1
#    else:
#        print("No match on: %s"%name)
#        nomatch+= 1
#
#   #print("%s %d"%(name,len(df1.loc[df1['StarName'] == name])))
#
#print("%d %d"%(match,nomatch))
#--------------

# Find the matching coordinates because the names don't always match
coord2 = SkyCoord(moretable['RA_d'],moretable['DEC_d'])
coord1 = SkyCoord(tconcat['RA2000'],tconcat['DEC2000'],unit=(u.hour,u.deg))
idx,d2d,d3d = coord1.match_to_catalog_sky(coord2)

#print(len(coord1),len(coord2))
for i in range(len(idx)):
    z = idx[i]
    #print("%s %s %s %.3f %.3f %%.3f"%(tconcat['StarName'][i],moretable['StarName'][z],coord1[i].to_string(),moretable['RA_d'][z],moretable['DEC_d'][z]),d2d[i])
    moretable['StarName'][z] = tconcat['StarName'][i]

tfinal = join(tconcat,moretable)

#--------------
# cute but not necessary.
# Make a coordinate pair from columns that look like coordinates :-)
#See http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord.guess_from_table
# an add a new column to table.
#tconcat['coordinates'] = SkyCoord.guess_from_table(t,unit=(u.hourangle,u.degree))
#--------------


# Astropy IPAC format does not allow non-alphanumeric characters in column names, even with non-strict checking (DBMS=False),
# even though the actual IPAC standard allows them! http://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html
formats = dict()
for c in tfinal.colnames:
    s = c.replace("'","").replace("-","_").replace(" ","_")
    tfinal[c].name = s
    # force a sensible output format for floats
    if tfinal[s].dtype == np.dtype('float64'): formats[s] = "%0.3f"
tfinal.write("sdss_standards.tab",format='ipac',DBMS=False,formats=formats,overwrite=True)
tfinal.write("sdss_standards.votab",format='votable',overwrite=True)

# Now make a SEDFitter source input file, in its peculiar format

cols = ["StarName","RA_d","DEC_d","u","rms_u","g","rms_g","r","rms_r","i","rms_i","z","rms_z",
        "FLUX_U","FLUX_ERROR_U","FLUX_B","FLUX_ERROR_B","FLUX_V","FLUX_ERROR_V",
        "FLUX_R","FLUX_ERROR_R","FLUX_I","FLUX_ERROR_I",
        "FLUX_J","FLUX_ERROR_J","FLUX_H","FLUX_ERROR_H","FLUX_K","FLUX_ERROR_K"
       ]

fsm = fm.FilterSetManager()
for i in tfinal:
