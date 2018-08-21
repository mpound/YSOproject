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
do_show = False
do_fit  = True
fudge_Gaia = False
GAIA_FUDGE_FACTOR = 5.0
fudge_Spitzer= False
SPITZER_FUDGE_FACTOR = 1./5.0

infile = "StarsWithAllPhotometry.votab"
tsources = Table.read(infile)

# Fix gaia errors lower limit is 0.004 (See GAIADR2 paper)
for b in ['phot_g_mean_mag_error', 'phot_bp_mean_mag_error', 'phot_rp_mean_mag_error']:
    tsources[b] = np.clip(tsources[b],a_min=0.004,a_max = None)
#tsources['phot_g_mean_mag_error']  = np.clip(tsources['phot_g_mean_mag_error'],a_min=0.004,a_max=None)
#tsources['phot_bp_mean_mag_error']  = np.clip(tsources['phot_bp_mean_mag_error'],a_min=0.004,a_max=None)
#tsources['phot_rp_mean_mag_error']  = np.clip(tsources['phot_rp_mean_mag_error'],a_min=0.004,a_max=None)

# clip SDSS also
for b in ["u_SDSS_err", "g_SDSS_err", "r_SDSS_err","i_SDSS_err", "z_SDSS_err"]:
    tsources[b] = np.clip(tsources[b],a_min=0.01,a_max = None)

if fudge_Gaia:
    for b in ['phot_g_mean_mag_error', 'phot_bp_mean_mag_error', 'phot_rp_mean_mag_error']:
       tsources[b]= tsources[b]*GAIA_FUDGE_FACTOR
if fudge_Spitzer:
    for b in [ "I1_flux_err", "I2_flux_err", "I3_flux_err" , "I4_flux_err"  ,      "h_msigcom","ks_msigcom"]:
       tsources[b] = tsources[b]*SPITZER_FUDGE_FACTOR

photometry_cols = [
        (fm.SDSS_u,"u_SDSS","u_SDSS_err"),
        (fm.SDSS_g,"g_SDSS","g_SDSS_err"),
#        (fm.SDSS_r,"r_SDSS","r_SDSS_err"),
#        (fm.SDSS_i,"i_SDSS","i_SDSS_err"),
#        (fm.SDSS_z,"z_SDSS","z_SDSS_err"),
        (fm.GAIA_G2,"phot_g_mean_mag","phot_g_mean_mag_error"),
        (fm.GAIA_B2,"phot_bp_mean_mag","phot_bp_mean_mag_error"),
        (fm.GAIA_R2,"phot_rp_mean_mag","phot_rp_mean_mag_error"),
        (fm.TWOMASS_J,"j_m","j_msigcom"),
        (fm.TWOMASS_H,"h_m","h_msigcom"),
        (fm.TWOMASS_K,"ks_m","ks_msigcom"),
        (fm.IRAC1,"I1_flux","I1_flux_err"),
        (fm.IRAC2,"I2_flux","I2_flux_err"),
        (fm.IRAC3,"I3_flux","I3_flux_err"),
        (fm.IRAC4,"I4_flux","I4_flux_err"),
        (fm.MIPS1,"M1_flux","M1_flux_err"),
#        (fm.MIPS2,"M2_flux","M2_flux_err")
       ]
distance_cols = ['Distance','Distance_err']
av_cols       = ['GB_Av','GB_Av_err']

outdir="fits_GAIA_stars_extinction_noriz/"
plotdir="plots_GAIA_stars_extinction_noriz/"
badfile = 'badfits.GAIA_stars_extinction'
#outdir="fits_GAIA_stars_extinction__fudge_gaia_noriz/"
#plotdir="plots_GAIA_stars_extinction_fudge_gaia_noriz/"
if fudge_Gaia and fudge_Spitzer:
    outdir="fits_GAIA_stars_extinction__fudge_gaia_spitzer_noriz_5/"
    plotdir="plots_GAIA_stars_extinction_fudge_gaia_spitzer_noriz_5/"

#outdir="fits_GAIA_stars_extinction_no_mips_fudge_gaia/"
#plotdir="plots_GAIA_stars_extinction_no_mips_fudge_gaia/"
#outdir="fits_GAIA_stars_extinction_fudge_gaia/"
#plotdir="plots_GAIA_stars_extinction_fudge_gaia/"
#skipme = ["M1_flux"]
skipme = [""]

photplotdir="plots_GAIA_stars_photometry/"
for z in [outdir, plotdir,photplotdir]:
    if not os.path.exists(z):
            os.mkdir(z)

#outdir='fits_ugrizJHK_extinction/'
#plotdir="plots_ugrizJHK_extinction/"
#outdir='fits_ugrizUBVRIJHK/'
#plotdir="plots_ugrizUBVRIJHK/"
#outdir='fits_UBVRIJHK_extinction/'
#plotdir="plots_UBVRIJHK_extinction/"
#outdir='fits_UBVRI/'
#plotdir="plots_UBVRI/"
#outdir='fits_ugriz/'
#plotdir="plots_ugriz/"
#outdir = "fits_ugriz_extinction/"
#plotdir = "plots_ugriz_extinction/"
#outdir = "fits_gaia2jhk_extinction/"
#plotdir = "plots_gaia2jhk_extinction/"

#shutil.rmtree(fitdir)
#shutil.rmtree(plotdir)

apertures = []
filters   = []
fsm = fm.FilterSetManager()
for i in range(len(photometry_cols)):
    apertures.append(3.0)
    #tel = fm._valid_bands[cols[i][0]].lower()
    #wave=fsm.wavelength(tel,cols[i][0]).to(u.micron)
    filters.append(photometry_cols[i][0])
    #filters.append(wave.value)

#MIPS Apertures
#apertures[-1] = 18.0
apertures[-1] = 6.0
apertures *= u.arcsec
#filters   *= u.micron
print(filters)
print(apertures)
colors = [
 #          u            g          r      i       z        
#           "violet","lightgreen","magenta","indianred","cyan", 
           "violet","violet","violet","violet","violet",
#           "black","black","black", #"black","black",
 #          G      GB      GR
            
          #"darkgreen","turquoise","plum",
          "darkgreen","darkgreen","darkgreen",
#           "black","black","black",
#            U     B       V       R       I        
#"black","blue","green","red","indigo",
#           "black","black","black","black","black",
#            J    H     K
#          "orange","pink","darkred"
#           "black","black","black",
           "blue","blue","blue",
# IRAC/MIPS
           "red","red","red","red","red"
         ]
markers = [
 #          ugriz        
          "+","+","+","+","+", 
 #          G      GB      GR
            
          "o","o","o",
#            U     B       V       R       I        
#          "v","v","v","v","v",
#            J    H     K
          "s","s","s",
# SPITZER I1 I2 I3 I4     M1
          "^","^","^","^",
          "*"
          ]
lu = dict()
lu['o'] = 'GAIA'
lu['+'] = 'ugriz'
#lu['v'] = 'UBVRI'
lu['^'] = 'IRAC'
lu['s'] = 'JHK'
lu['*'] = 'MIPS'

seds = []
for i in tsources:
    abad = False
    d = i[distance_cols[0]]
    de = i[distance_cols[1]]
#    print("ACTUAL DISTANCE = %s +/- %s"%(d*u.pc,de*u.pc))
    av  = i[av_cols[0]]
    ave  = i[av_cols[1]]
    if ma.is_masked(d) or ma.is_masked(de) :
            continue

    #print(type(av),av)
    s = SED(i['source_id'],d*u.pc,de*u.pc,i['ra']*u.degree,i['dec']*u.degree,av=av,averr=ave)
    for c in photometry_cols:
        if c[1] in skipme: validity = 0
        else: validity=1
        if np.ma.is_masked(i[c[1]]) or np.ma.is_masked(i[c[2]]):
#            print("flux and/or error is masked for Source %s band %s...setting validity to zero."%(s._name,c[0]))
#            validity = 0
            print("flux and/or error is masked for Source %s band %s...skipping ."%(s._name,c[0]))
            abad = True
        else:
            if tsources[c[1]].unit=='mag':
                s.addData(c[0],u.Magnitude(i[c[1]]),u.Magnitude(i[c[2]]),validity)
            else:
                s.addData(c[0],i[c[1]]*u.mJy,i[c[2]]*u.mJy,validity)
    s.set_upper_limits()
    if not abad: 
        seds.append(s)
        if do_plot:
            for x,y,c,m in zip(s.wavelengths(),s.fluxes(),colors,markers):
                plt.scatter(x,y,c=c,marker=m,s=20)
            #for x,y,m in zip(s.wavelengths(),s.fluxes(),markers):
                #plt.scatter(x,y,c='k',marker=m,s=20)
            #print(s.wavelengths(),s.fluxes())
            mymin = np.min(s.fluxes())
            mymax = np.max(s.fluxes())
            plt.gca().set_ylim(0.1*mymin,10*mymax)
            plt.title(s._name)
            patches = [plt.plot([],[],marker=z,ls="",color='k',label=lu[z])[0] for z in lu.keys()]
            plt.legend(handles=patches,loc='upper right',numpoints=1,ncol=1, framealpha=0.25) #bbox_to_anchor=(.9,0.65),framealpha=0.25)
            plt.gca().set_xscale('log')
            plt.gca().set_yscale('log')
            plt.savefig(photplotdir+s._name+"_Photometry.png")
            if do_show: plt.show()

#print(s.sedfitterinput())
print("Starting with %d sources" % len(seds))


# ugh sedfitter.fit() only allows one distance so we have to do 
# each source one at a time.
default_av_range=[0., 5.]
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


bad = open(badfile,"w")
if do_fit:
    for s in seds:
        #source = Source.from_ascii(s.sedfitterinput())
        thisbad=False
        datafile='tmpdir/'+s._name+"_input.txt"
        f = open(datafile,"w")
        f.write(s.sedfitterinput())
        f.close()
        #print(s.sedfitterinput())
        #nospacename = s._name.replace(' ','_')
        nospacename = str(s._name).replace(' ','_')
        outfit = outdir+"/"+nospacename+'.sedfit'
        distance_range = u.pc*s.distrange_pc()
        print("DISTANCE RANGE is ",distance_range)
        if ma.is_masked(s._av):
           av_range = default_av_range
        else:
           av_range = s.avrange()
        if None in av_range:
            av_range = default_av_range
        #av_range = default_av_range
        print("AV RANGE is ",av_range)
        select_format=('N',25)
        try:
            fit(datafile,filter_names=filters,
                    apertures=apertures,extinction_law=extinction,
                    distance_range=distance_range,
                    av_range=av_range,model_dir=sed_model_dir,
                    output_convolved=True, output=outfit,
                    remove_resolved=False,output_format=select_format)
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
            write_parameters(outfit,outdir+"/"+nospacename+"_parameters.txt",select_format=select_format)
            write_parameter_ranges(outfit,outdir+"/"+nospacename+"_parameters_ranges.txt",select_format=select_format)
            plot(input_fits=outfit, output_dir=plotdir+nospacename,format='png',
                     plot_mode='A',plot_name=True,
                     select_format=select_format,
                     show_convolved=False, show_sed=True)

bad.close()
