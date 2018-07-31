#!/usr/bin/env python
# REQUIRES Python 3
#
# Marc Pound - mpound@umd.edu - Thu Jul 12 15:39:34 EDT 2018
#
# Test fit of new filters.  This follows the example at
# http://sedfitter.readthedocs.io/en/stable/fitting.html#fitting-the-models-to-the-data
#
# If you get a  
#ValueError: operands could not be broadcast togethe
#make sure number of filters matches number of points
# however, there must be convolved models for every filter even if unused

import os
import glob
import sys
import traceback
import numpy as np
import filtermanage as fm

from astropy import units as u

from sedfitter import (fit, plot, plot_params_1d, plot_params_2d,
                       write_parameters, write_parameter_ranges)
from sedfitter.source import Source
from sedfitter.filter import Filter
from sedfitter.extinction import Extinction

# Define path to models
topdir        = '/n/subaruraid/mpound/'
model_dir     = topdir+'sedfittermodels/'
sed_model_dir = model_dir+'models_r17/s---smi'
dust_model    = model_dir+'whitney.r550.par'
data_dir      = topdir+'SFM/'
#datafile   = data_dir+'source_sedfit_input'
#datafile      = data_dir+'source_input.use'
onedatafile      = data_dir+'example_one_source.joncour'
fulldatafile      = data_dir+'input_sources.joncour'
shortdatafile      = data_dir+'shortsources.joncour'
visibledatafile    = data_dir+'example_visible_source.joncour'
distance      = 436.0
disterr       = 9.0
distance_range = [distance-disterr,distance+disterr] * u.pc
av_range=[0., 40.]


# Read in extinction law. We read in columns 1 (the wavelength in microns) and
# 4 (the opacity in cm^2/g)
extinction = Extinction.from_file(dust_model, columns=[0, 3],
                                  wav_unit=u.micron, chi_unit=u.cm**2 / u.g)

# Define filters and apertures
filters   = ['2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4', 'M1', 'M2', 'PACS1', 'PACS2']
apertures = [ 3.,    3.,   3.,   3.,   3.,   3.,   3.,  3.,    3., 3.,     3.] * u.arcsec

filters   = ['2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4']
apertures = [ 3.,    3.,   3.,   3.,   3.,   3.,   3.]*u.arcsec

if False:
    filters   = ['2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4', 'G','BP','RP','u','g','r','i','z']
    apertures = [ 3.,    3.,   3.,   3.,   3.,   3.,   3., 3., 3., 3., 3., 3., 3., 3., 3.] * u.arcsec
    ext = "_all_filts"
else:
    ext = ""
#filters   = ['2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4', 6735.4*u.angstrom, 5319.9*u.angstrom, 7992.9*u.angstrom, 3561.8*u.angstrom, 4718.9*u.angstrom, 6185.2*u.angstrom, 7499.7*u.angstrom, 8961.5*u.angstrom]
#apertures = [ 3.,    3.,   3.,   3.,   3.,   3.,   3., 3., 3., 3., 3., 3., 3., 3., 3.] * u.arcsec
# Run the fitting
ext="onesource_testr17"
outfile       = data_dir+'one_source_sedfit_output.info'+ext
if True:
    try:
        fit(data=onedatafile,
            filter_names=filters, 
            apertures=apertures, 
            model_dir=sed_model_dir,
            output=outfile+ext,
            extinction_law=extinction,
            distance_range=distance_range,
            av_range=av_range,
            output_convolved=False, remove_resolved=True)
    except ValueError, v:
        traceback.print_exc()
        print("\n### ERROR: %s"%v)
        print("### Check that number of filters matches number of photometric data points")
        sys.exit(255)


if False:
    mysources = dict()
    with open(visibledatafile) as srcfile:
        for line in srcfile:
            s = Source.from_ascii(line)
            mysources[s.name] = s
#    print mysources
#    myfilternames   = ['2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4', 'G','BP','RP','u','g','r','i','z']
#    waves = np.array([12350.0*u.angstrom, 16620.0*u.angstrom, 21590.0*u.angstrom, 35572.6*u.angstrom, 45049.3*u.angstrom, 57385.7*u.angstrom, 79273.7*u.angstrom, 6735.4*u.angstrom, 5319.9*u.angstrom, 7992.9*u.angstrom, 3561.8*u.angstrom, 4718.9*u.angstrom, 6185.2*u.angstrom, 7499.7*u.angstrom, 8961.5*u.angstrom],dtype=type(1.0*u.angstrom))
    myfilternames   = ['GAIAG','GAIAB','GAIAR','SLOANu','SLOANg','SLOANr','SLOANi','SLOANz', 'M1']
    waves = np.array([6735.4*u.angstrom, 5319.9*u.angstrom, 7992.9*u.angstrom, 3561.8*u.angstrom, 4718.9*u.angstrom, 6185.2*u.angstrom, 7499.7*u.angstrom, 8961.5*u.angstrom, 238433.1*u.angstrom],dtype=type(1.0*u.angstrom))
#    nus = waves.to(u.Hz,equivalencies=u.spectral())
#    resp = np.ndarray([1])
#    resp[0] = 1.0
    myfilters = []
    for i in range(len(myfilternames)):
        thisfilter = dict()
#        nu1 =  np.array([waves[i].to(u.Hz,equivalencies=u.spectral())],dtype=type(waves[i]))
        #nu = np.array( waves.to(u.Hz,equivalencies=u.spectral()))
        #g = Filter(myfilternames[i], waves[i], None, resp)
        thisfilter['name']=myfilternames[i]
        thisfilter['wav']=waves[i]
        thisfilter['aperture_arcsec']=3.0
        myfilters.append(thisfilter)

#mywav = np.array([f['wav'].to(u.micron).value for f in filters])


# For the remaining commands, we always select the models with chi^2-chi_best^2
# per datapoint less than 3.
select_format = ('F', 10)
select_format = ('F', 3)

# Make SED plots

ext="onesource_nojhk+mips"

labels = []
labels.append("Best Fit: A$_{\\rm V}=9.7$, M$_{\\rm disk}=1.4\\times 10^{-4} {\\rm M}_{\\odot}$, $\\dot{{\\rm M}}=0$, T$_{\\rm eff}=4000$, L$= 0.5{\\rm L_{\odot}}$")
labels.append("w/Optical+M24 data: A$_{\\rm V}\\sim 2$, M$_{\\rm disk}\\sim 10^{-6} {\\rm M}_{\\odot}$, $\dot{{\\rm M}}\\sim 5\\times 10^{-8}$, T$_{\\rm eff}=4000$, L$= 0.25{\\rm L_{\\odot}}$")
outfile       = data_dir+'one_source_sedfit_output.info'+ext
plot(outfile, 'plots_seds'+ext, plot_max=100, select_format=select_format,mysources=mysources, myfilters=myfilters,y_mode='M',y_range=(6E-14,5E-10),x_mode='M',x_range=(2E-1,80),plot_info=False,plot_name=False,mylabels=labels)

ext="onesource+mips"

labels = []
labels.append("Best Fit: A$_{\\rm V}=0$, M$_{\\rm disk}=5\\times 10^{-3} {\\rm M}_{\\odot}$, $\\dot{{\\rm M}}=2\\times 10^{-7}$, T$_{\\rm eff}=4000$, L$= 3.0{\\rm L_{\\odot}}$")
labels.append("w/Optical+M24 data: A$_{\\rm V}\\sim 2$, M$_{\\rm disk}\\sim 10^{-6} {\\rm M}_{\\odot}$, $\\dot{{\\rm M}}\\sim 5\\times 10^{-8}$, T$_{\\rm eff}=4000$, L$= 0.25{\\rm L_{\\odot}}$")
outfile       = data_dir+'one_source_sedfit_output.info'+ext
plot(outfile, 'plots_seds'+ext, plot_max=100, select_format=select_format,mysources=mysources, myfilters=myfilters,y_mode='M',y_range=(6E-14,5E-10),x_mode='M',x_range=(2E-1,80),plot_info=False,plot_name=False,mylabels=labels)

if False:
    # Make histograms of the disk mass
    plot_params_1d(outfile, 'MDISK', data_dir+'plots_mdisk'+ext,
                   log_x=True, select_format=select_format)

    # Make 2-d plots of the envelope infall rate vs disk mass
    plot_params_2d(outfile, 'MDISK', 'MDOT', data_dir+'plots_mdot_mdisk'+ext,
                   log_x=True, log_y=True, select_format=select_format)

    # Write out all models with a delta chi^2-chi_best^2 per datapoint < 3
    write_parameters(outfile, data_dir+'parameters.txt'+ext,
                     select_format=select_format)

    # Write out the min/max ranges corresponding to the above file
    write_parameter_ranges(outfile, data_dir+'parameter_ranges.txt'+ext,
                           select_format=select_format)
