#!/usr/bin/env python
#
# create a SEDFitter input file from Serpens Cross-match data 
# provided by Isabelle
#
import numpy as np
import numpy.ma as ma
from astropy.io import ascii
from astropy.table import Table
import filtermanage as fm

#              JHK    I1-4    G BP RP ugriz  
allfilts = True
useIR    = True
if allfilts:
    nfilts =        3       +4      +3      +5
    use_filts = [ 1,1,1, 1,1,1,1, 1,1,1, 1,1,1,1,1 ] 
else:
    nfilts =        3       +4     
    use_filts = [ 1,1,1, 1,1,1,1]
#f = Table.read('/n/subaruraid/mpound/ADAP/XMatch_SerpMain_28sources.csv',format='ascii.csv')
f = Table.read('XMatch_SerpMain_28sources_updated.csv',format='ascii.csv')
fsm = fm.FilterSetManager()
for i in f:
    k = 0
    line = "%s  %f  %f" % ( i['Source_Get'],i['RAJ2000_Get'],i['DEJ2000_Get'])
    # JHK
    photline=""
    if useIR:
        photline+="%5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e " % (
            fsm.magtoflux(fm.TWOMASS,fm.TWOMASS_K,i['Jmag_Get']).value,
            fsm.magtoflux(fm.TWOMASS,fm.TWOMASS_K,i['Jmag_Get']).value*(i['e_Jmag_Get'] if i['e_Jmag_Get']>0 else 0.01), 
            fsm.magtoflux(fm.TWOMASS,fm.TWOMASS_H,i['Hmag_Get']).value,
            fsm.magtoflux(fm.TWOMASS,fm.TWOMASS_H,i['Hmag_Get']).value*(i['e_Hmag_Get'] if i['e_Hmag_Get']>0 else 0.01), 
            fsm.magtoflux(fm.TWOMASS,fm.TWOMASS_K,i['Ksmag_Get']).value,
            fsm.magtoflux(fm.TWOMASS,fm.TWOMASS_K,i['Ksmag_Get']).value*(i['e_Ksmag_Get'] if i['e_Ksmag_Get']>0 else 0.01),
        #SPITZER
            fsm.magtoflux(fm.SPITZER,fm.IRAC1,i['3.6mag_Get']).value,
            fsm.magtoflux(fm.SPITZER,fm.IRAC1,i['3.6mag_Get']).value*(i['e_3.6mag_Get'] if i['e_3.6mag_Get'] else 0.01), 
            fsm.magtoflux(fm.SPITZER,fm.IRAC2,i['4.5mag_Get']).value,
            fsm.magtoflux(fm.SPITZER,fm.IRAC2,i['4.5mag_Get']).value*(i['e_4.5mag_Get'] if i['e_4.5mag_Get'] else 0.01),
            fsm.magtoflux(fm.SPITZER,fm.IRAC3,i['5.8mag_Get']).value,
            fsm.magtoflux(fm.SPITZER,fm.IRAC3,i['5.8mag_Get']).value*(i['e_5.8mag_Get'] if i['e_5.8mag_Get'] else 0.01), 
            fsm.magtoflux(fm.SPITZER,fm.IRAC4,i['8.0mag_Get']).value,
            fsm.magtoflux(fm.SPITZER,fm.IRAC4,i['8.0mag_Get']).value*(i['e_8.0mag_Get'] if i['e_8.0mag_Get'] else 0.01))
        k=6
    else:
        k=-1
        use_filts = [1,1,1, 1,1,1,1,1]
        nfilts =        3       +5     

    if  allfilts:
        #fm.GAIA
        k+=1
        val = fsm.magtoflux(fm.GAIA,fm.GAIA_G,i['phot_g_mean_mag_Ga']).value
        if ma.is_masked(val): 
              photline += " -999  -999 "
              use_filts[k] = 0
        else: 
              photline+=" %5.4e  %5.4e " % (val, val*i['phot_g_mean_flux_error_Ga']/i['phot_g_mean_flux_Ga'])
              use_filts[k] = 1
        k+=1
        val = fsm.magtoflux(fm.GAIA,fm.GAIA_B,i['phot_bp_mean_mag_Ga']).value
        if ma.is_masked(val): 
              photline += " -999  -999 "
              use_filts[k] = 0
        else: 
              photline+=" %5.4e  %5.4e " % (val, val*i['phot_bp_mean_flux_error_Ga']/i['phot_bp_mean_flux_Ga'])
              use_filts[k] = 1
        k+=1
        val = fsm.magtoflux(fm.GAIA,fm.GAIA_R,i['phot_rp_mean_mag_Ga']).value
        if ma.is_masked(val): 
              photline += " -999  -999 "
              use_filts[k] = 0
        else: 
              photline+=" %5.4e  %5.4e " % (val, val*i['phot_rp_mean_flux_error_Ga']/i['phot_rp_mean_flux_Ga'])
              use_filts[k] = 1

        #SLOAN
        k+=1
        val = fsm.magtoflux(fm.SDSS,fm.SDSS_u,i['u_mag_Ga']).value
        if ma.is_masked(val):
              photline += " -999  -999 "
              use_filts[k] = 0
        else: 
              photline+=" %5.4e  %5.4e " % (val, val*i['u_mag_error_Ga'])
              use_filts[k] = 1
        k+=1
        val = fsm.magtoflux(fm.SDSS,fm.SDSS_g,i['g_mag_Ga']).value
        if ma.is_masked(val):
              photline += " -999  -999 "
              use_filts[k] = 0
        else: 
              photline+=" %5.4e  %5.4e " % (val, val*i['g_mag_error_Ga'])
              use_filts[k] = 1
        k+=1
        val = fsm.magtoflux(fm.SDSS,fm.SDSS_r,i['r_mag_Ga']).value
        if ma.is_masked(val):
              photline += " -999  -999 "
              use_filts[k] = 0
        else: 
              photline+=" %5.4e  %5.4e " % (val, val*i['r_mag_error_Ga'])
              use_filts[k] = 1
        k+=1
        val = fsm.magtoflux(fm.SDSS,fm.SDSS_i,i['i_mag_Ga']).value
        if ma.is_masked(val):
              photline += " -999  -999 "
              use_filts[k] = 0
        else: 
              photline+=" %5.4e  %5.4e " % (val, val*i['i_mag_error_Ga'])
              use_filts[k] = 1
        k+=1
        val = fsm.magtoflux(fm.SDSS,fm.SDSS_z,i['z_mag_Ga']).value
        if ma.is_masked(val): 
              photline += " -999  -999 "
              use_filts[k] = 0
        else: 
              photline+=" %5.4e  %5.4e " % (val, val*i['z_mag_error_Ga'])
              use_filts[k] = 1

    filtsline = ""
    for j in range(nfilts): filtsline += " %d "%use_filts[j]
    line += filtsline
    line += photline

    print(line)
