#!/usr/bin/env python
import numpy as np
import numpy.ma as ma
from astropy.io import ascii
from astropy.table import Table,vstack
from astropy import units as u
from astropy.visualization import quantity_support
quantity_support()
import matplotlib.pyplot as plt
import pandas as pd

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

# you access a value by it's column name, so e.g. the first source name:
#print(t['c2d_ID'][0])
#mask1 = t['Prob_Galc'] < -0.6 
#mask2 = t['Prob_Galc'] > -1.47
#mask3 = np.logical_and(mask1,mask2) 
## this would have both masks applied.
# something not right because there are only 67 sources after this operation
#t[mask3]

GBysodir='/n/subaruraid/lgm/Spitzer-Serpens/GB-YSOc/'
GBfulldir='/n/subaruraid/lgm/Spitzer-Serpens/GB-Full/'
ysofname = 'catalog-SERAQU_%d-YSOc.tbl'
fullfname = 'catalog-SERAQU_%d-FULL.tbl'
# There are 6 regions. YSOc Region 5 is the area near Serpens East that we 
# have selected the one square degree from.  Similary Gb-Full Region 6.
indx = range(1,7)

plt.xlim(280.5,278.5)
plt.ylim(-0.5,1)
plt.axes().set_aspect('equal')

concat = False
if concat:
    ysotables = list()
    for i in indx:
        print(GBysodir+ysofname%(i))
        f=GBysodir+ysofname%(i)
        t =Table.read(f,format='ipac')
        ysotables.append(t)
        #plt.plot(t['ra']*u.degree,t['dec']*u.degree)

    # concatenate all the tables
    allyso = vstack(ysotables)
    plt.scatter(ysotables[5]['ra'],ysotables[5]['dec'],color='orange',s=3)
    df = allyso.to_pandas()
    df = df[(df.ra <= 280.25)& (df.ra >= 279.25) & (df.dec >= -0.25) & (df.dec <= 0.75)]
    mytable = Table.from_pandas(df)
    print("Serpens East YSOc sources: %d"%len(mytable))
    mytable.write("catalog-SerpensEast_GB_YSOc.tbl",format="ipac")
    df = df[(df.Prob_Galc >= -0.6)]
    print("Selected serpens east %d lines  PG >= -0.6" % len(df))
else:
    yso = Table.read("catalog-SerpensEast_GB_YSOc.tbl",format="ipac")
    print("Serpens East YSOc sources: %d {orange}"%len(yso))
    df = yso.to_pandas()
    df = df[(df.Prob_Galc >= -0.6)]
    print("Selected serpens east %d lines  PG >= -0.6" % len(df))
    plt.scatter(yso['ra'],yso['dec'],color='orange',s=3)


readbigtable = False
if readbigtable:
    #610K lines!
    GBSerpEast = GBfulldir+fullfname%(6)
    bigtable = Table.read(GBSerpEast,format='ipac')
    print("Read %d lines" % len(bigtable))
    # make a pandas dataframe to more easily select the data we want
    df=bigtable.to_pandas()
    xdf = df[(df.Prob_Galc > -1.47) & (df.Prob_Galc < -0.6)]
    print("Selected %d lines [-1.47 < PG < -0.6" % len(xdf))
    smalltable = Table.from_pandas(xdf)
    ascii.write(smalltable,'catalog-SERAQU_6_ProbGalc_-1.47to-0.6.tbl',format='ipac')
    xxdf = df[(df.Prob_Galc >= -0.6)]
    xsmalltable = Table.from_pandas(xxdf)
    print("Selected %d lines [PG >= -0.6]" % len(xxdf))
    ascii.write(xsmalltable,'catalog-SERAQU_6_ProbGalc_ge_-0.6.tbl',format='ipac')
else:
    xsmalltable = Table.read('catalog-SERAQU_6_ProbGalc_ge_-0.6.tbl',format='ipac')
    smalltable = Table.read('catalog-SERAQU_6_ProbGalc_-1.47to-0.6.tbl',format='ipac')
    print("Selected %d lines [-1.47 < PG < -0.6] {blue}" % len(smalltable))
    print("Selected %d lines [PG >= -0.6] {red - galaxies!}" % len(xsmalltable))
plt.scatter(smalltable['ra'],smalltable['dec'],color='blue',s=25,facecolor='none')
plt.scatter(xsmalltable['ra'],xsmalltable['dec'],color='red',s=25,facecolor='none')
plt.show()

