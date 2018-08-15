#!/usr/bin/env python
import os
import glob
import sys
import shutil
import traceback
import numpy as np
import numpy.ma as ma
import filtermanage as fm

from astropy import units as u

from sedfitter import (fit, plot, plot_params_1d, plot_params_2d,
                       write_parameters, write_parameter_ranges, Fitter, FitInfoFile)
from sedfitter.source import Source
from sedfitter.filter import Filter
from sedfitter.extinction import Extinction
from sedfitter.sed import SEDCube
from astropy.visualization import quantity_support
quantity_support()
import matplotlib.pyplot as plt

from sed import SED
import filtermanage as fm
from astropy.io import ascii
from astropy.table import Table


do_plot = False
do_fit  = True

tsources = Table.read("sdss_standards+gaia.votab")

cols = [
#        (fm.SDSS_u,"u","rms_u"),
#        (fm.SDSS_g,"g","rms_g"),
#        (fm.SDSS_r,"r","rms_r"),
#        (fm.SDSS_i,"i","rms_i"),
#        (fm.SDSS_z,"z","rms_z"),
        (fm.GAIA_G2,"phot_g_mean_mag","phot_g_mean_mag_error"),
        (fm.GAIA_B2,"phot_bp_mean_mag","phot_bp_mean_mag_error"),
        (fm.GAIA_R2,"phot_rp_mean_mag","phot_rp_mean_mag_error"),
#        (fm.BESSELL_U,"FLUX_U","FLUX_ERROR_U"),
#        (fm.BESSELL_B,"FLUX_B","FLUX_ERROR_B"),
#        (fm.BESSELL_V,"FLUX_V","FLUX_ERROR_V"),
#        (fm.BESSELL_R,"FLUX_R","FLUX_ERROR_R"),
#        (fm.BESSELL_I,"FLUX_I","FLUX_ERROR_I"),
        (fm.TWOMASS_J,"FLUX_J","FLUX_ERROR_J"),
        (fm.TWOMASS_H,"FLUX_H","FLUX_ERROR_H"),
        (fm.TWOMASS_K,"FLUX_K","FLUX_ERROR_K")
       ]
apertures = []
filters   = []
fsm = fm.FilterSetManager()
for i in range(len(cols)):
    apertures.append(5.0)
    #tel = fm._valid_bands[cols[i][0]].lower()
    #wave=fsm.wavelength(tel,cols[i][0]).to(u.micron)
    filters.append(cols[i][0])
    #filters.append(wave.value)
apertures *= u.arcsec
#filters   *= u.micron
print(filters)
print(apertures)
colors = [
 #          u            g          r      i       z        
#           "violet","lightgreen","magenta","indianred","cyan", 
#           "black","black","black","black","black",
 #          G      GB      GR
            
#          "darkgreen","turquoise","plum",
           "black","black","black",
#            U     B       V       R       I        
#          "black","blue","green","red","indigo",
#           "black","black","black","black","black",
#            J    H     K
#          "orange","pink","darkred"
           "black","black","black",
         ]
markers = [
 #          ugriz        
#           "+","+","+","+","+", 
 #          G      GB      GR
            
          "o","o","o",
#            U     B       V       R       I        
#          "v","v","v","v","v",
#            J    H     K
          "s","s","s"
          ]
lu = dict()
lu['+'] = 'ugriz'
lu['o'] = 'GAIA'
lu['v'] = 'UBVRI'
lu['s'] = 'JHK'
#outdir='fits_ugrizJHK/'
#plotdir="plots_ugrizJHK/"
#outdir='fits_ugrizUBVRIJHK/'
#plotdir="plots_ugrizUBVRIJHK/"
#outdir='fits_UBVRIJHK/'
#plotdir="plots_UBVRIJHK/"
#outdir='fits_UBVRI/'
#plotdir="plots_UBVRI/"
#outdir='fits_ugriz/'
#plotdir="plots_ugriz/"
#outdir = "fits_ugriz_extinction/"
#plotdir = "plots_ugriz_extinction/"
outdir = "fits_gaia2jhk_extinction/"
plotdir = "plots_gaia2jhk_extinction/"

seds = []
for i in tsources:
    bad = False
    if ma.is_masked(i['Distance_distance']) or ma.is_masked(i['Distance_merr']) or ma.is_masked(i['Distance_perr']): 
            continue

    av = 3.1*i['EB_V']
    #print(type(av),av)
    s = SED(i['StarName'],i['Distance_distance']*u.pc,(i['Distance_merr'],i['Distance_perr'])*u.pc,i['RA_d']*u.degree,i['DEC_d']*u.degree,av=av)
    for c in cols:
        validity=1
        if np.ma.is_masked(i[c[1]]) or np.ma.is_masked(i[c[2]]):
#            print("flux and/or error is masked for Source %s band %s...setting validity to zero."%(s._name,c[0]))
#            validity = 0
            print("flux and/or error is masked for Source %s band %s...skipping ."%(s._name,c[0]))
            bad = True
        else:
            s.addData(c[0],u.Magnitude(i[c[1]]),u.Magnitude(i[c[2]]),validity)
    if not bad: 
        seds.append(s)
        if do_plot:
            for x,y,c,m in zip(s.wavelengths(),s.fluxes(),colors,markers):
                plt.scatter(x,y,c=c,marker=m,s=20)
            plt.title(s._name)
            patches = [plt.plot([],[],marker=z,ls="",color='k',label=lu[z])[0] for z in lu.keys()]
            plt.legend(handles=patches,loc='best',numpoints=1,ncol=1, bbox_to_anchor=(.9,0.65),framealpha=1)
            plt.show()

#print(s.sedfitterinput())
print("Starting with %d sources" % len(seds))


# ugh sedfitter.fit() only allows one distance so we have to do 
# each source one at a time.
default_av_range=[0., 4.]
dust_model = 'whitney.r550.par'
topdir        = '/n/subaruraid/mpound/'
model_dir     = topdir+'sedfittermodels/'
#model_dir     = '/lupus3/mpound/filter_convolve/'
sed_model_dir = model_dir+'models_r17/s---s-i/'
extinction = Extinction.from_file(dust_model, columns=[0, 3],wav_unit=u.micron, chi_unit=u.cm**2 / u.g)

if False:
    sedmodels = SEDCube.read(sed_model_dir+"flux.fits")
    print(sedmodels.names)
    print(sedmodels.val.shape)
    s = sedmodels.get_sed('03ZZRVTe_01')
    plt.loglog(s.wav, s.flux.transpose(), 'k-', alpha=0.5)
    plt.loglog(seds[0].wavelengths(),seds[0].fluxes())
    plt.show()
#plt.ylim(1e-2, 1e8)


bad = open("badfits_gaiajhk","w")
if do_fit:
    for s in seds:
        #source = Source.from_ascii(s.sedfitterinput())
        thisbad=False
        datafile="tempsource.txt"
        f = open(datafile,"w")
        f.write(s.sedfitterinput())
        f.close()
        #print(s.sedfitterinput())
        nospacename = s._name.replace(' ','_')
        outfit = outdir+nospacename+'.sedfit'
        distance_range = [s._distance.value+s._disterr[0].value,s._distance.value+s._disterr[1].value]*u.pc
        if ma.is_masked(s._av):
           av_range = default_av_range
        else:
           if av<=0:
               av_range = [0 , av+1.5]
           else:
               av_range = [0 , av*2.0]
#        distance_range = [s._distance.value-50,s._distance.value+50]*u.pc
        #fitter = Fitter(filter_names=filters,
        try:
            fit(datafile,filter_names=filters,
                    apertures=apertures,extinction_law=extinction,
                    distance_range=distance_range,
                    av_range=av_range,model_dir=sed_model_dir,
                    output_convolved=False, output=outfit,
                    remove_resolved=False)
        except Exception as e:
            bad.write("%s  %s\n"% (nospacename,e))
            thisbad = True
            print(e)
        #info = fitter.fit(source)
        #FitInfo.write(info,nospacename+'.sedfit')
        if thisbad:
            continue
        else:
            #z = FitInfoFile(outdir+nospacename+'.sedfit',mode="r")
            ##print(type(z))
            #good.append(z)
            write_parameters(outfit,outdir+nospacename+"_parameters.txt",select_format=('F',10))
            plot(input_fits=outfit, output_dir=plotdir+nospacename,format='png',
                     plot_mode='A',plot_name=True,
                     select_format=('F',10),
                     show_convolved=False, show_sed=True)

bad.close()
