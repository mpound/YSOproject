#!/usr/bin/env python
import numpy as np
import numpy.ma as ma
from astropy.io import ascii
from astropy.table import Table,vstack
from astropy import units as u
from astropy.visualization import quantity_support
quantity_support()
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Read in the c2d Serpens full catalog.
# See catalog-SER-FULL.dd for the Data Dictionary that
# indicates what each column means

# The table is 377K lines long, this will take a while!
# This is c2d Serpens, we want gould's belt
#t = Table.read('catalog-SER-FULL.tbl',format='ipac')
# The available column names
#print(t.colnames)

# Some columns of interest
# c2d_ID       Source name 
#  ra          Right Ascension [J2000] 
# D_ra         Uncertainty in Right Ascension
# dec          Declination [J2000]  
# D_dec        Uncertainty in Declination 
# Prob_Galc    Non-normalized probability source is a galaxy 
#              YSO candidates have Prob_Galc > -0.6 (Orange in figure)
#              Blue source in figure have -1.47 < Prob_Galc < -0.6

datadir = 'GB-YSOc/'
ysofname = datadir+'catalog-SERAQU_6-YSOc.tbl'
ysotable = Table.read(ysofname,format='ipac')
ax = plt.gca()
plt.xlim(280.5,279.)
plt.ylim(-0.6,0.9)
plt.axes().set_aspect('equal')
plt.scatter(ysotable['ra'],ysotable['dec'],color='orange',s=3)

# Now sources with Prob_Galc > -1.47) AND (Prob_Galc < -0.6) 
secondarysources = datadir+'catalog-SERAQU_6_ProbGalc.tbl'
smalltable = Table.read(secondarysources,format='ipac')
plt.scatter(smalltable['ra'],smalltable['dec'],color='blue',s=3)

SelectData = False
if SelectData:
    # this selects the 1 square degree test patch
    df = smalltable.to_pandas()
    df = df[(df.ra <= 280.25)& (df.ra >= 279.25) & (df.dec >= -0.25) & (df.dec <= 0.75)]
    selectedtable = Table.from_pandas(df)
    selectedtable.write(datadir+"serpens_test_area.tbl",format="ipac")
else:
    selectedtable = Table.read(datadir+"serpens_test_area.tbl",format="ipac")

plt.scatter(selectedtable['ra'],selectedtable['dec'],color='green',s=25, facecolors='none')
# The 1 square degree box defined for this test
rect = patches.Rectangle((280.25,-0.25),-1.0,1.0,fill=False,color='green')
ax.add_patch(rect)
plt.show()

