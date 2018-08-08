#!/usr/bin/env python
from sed import SED
import filtermanage as fm
from astropy.io import ascii
from astropy.table import Table,join,Column
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import quantity_support
quantity_support()  
import matplotlib.pyplot as plt

cols = [
        (fm.SDSS_u,"u","rms_u"),
        (fm.SDSS_g,"g","rms_g"),
        (fm.SDSS_r,"r","rms_r"),
        (fm.SDSS_i,"i","rms_i"),
        (fm.SDSS_z,"z","rms_z"),
#        (fm.BESSEL_U,"FLUX_U","FLUX_ERROR_U"),
#        (fm.BESSEL_B,"FLUX_B","FLUX_ERROR_B"),
#        (fm.BESSEL_V,"FLUX_V","FLUX_ERROR_V"),
#        (fm.BESSEL_R,"FLUX_R","FLUX_ERROR_R"),
#        (fm.BESSEL_I,"FLUX_I","FLUX_ERROR_I"),
#        (fm.TWOMASS_J,"FLUX_J","FLUX_ERROR_J"),
#        (fm.TWOMASS_H,"FLUX_H","FLUX_ERROR_H"),
#        (fm.TWOMASS_K,"FLUX_K","FLUX_ERROR_K")
       ]

#print(cols)

tfinal = Table.read("sdss_standards.votab")
seds = []
for i in tfinal:
    s = SED(i['StarName'],i['Distance_distance']*u.pc,(i['Distance_merr'],i['Distance_perr'])*u.pc,i['RA_d']*u.degree,i['DEC_d']*u.degree)
    for c in cols:
        s.addData(c[0],u.Magnitude(i[c[1]]),u.Magnitude(i[c[2]]),1)
    seds.append(s)
print(seds[0].header())
counter = 0
#spec = gridspec.Gridspec(nrows=4,ncols=4)
if False:
    for s in seds:
         if counter % 16 == 0:
            ax1=ax2=0
            fig,ax = plt.subplots(ncols=4,nrows=4)
         axs=ax.reshape(-1)
    #    print(s.sedfitterinput())
         axs[counter].scatter(s.wavelengths(),s.fluxes(),c='k-')
         axs[counter].set_title(s._name)
         if counter % 16 == 15:
            plt.show()
            counter = -1
         counter += 1
         
         
for s in seds:
    plt.scatter(s.wavelengths(),s.fluxes())
    plt.title(s._name,pad=-12)
    plt.show()



