#!/usr/bin/env python
import numpy as np
from astropy.io import ascii
from astropy.table import Table, Column, vstack
import astropy.units as u
from astropy.constants import c
from sedfitter.sed import SEDCube
from sedfitter.convolved_fluxes import ConvolvedFluxes
from astropy.visualization import quantity_support
quantity_support()
import matplotlib.pyplot as plt

class AlphaMaker:
    def __init__(self,model_dir=None,model_type=None):
        if model_dir == None:
           self._model_dir = "/n/subaruraid/mpound/sedfittermodels/models_r17/"
        else:
           self._model_dir = model_dir
        if model_type == None:
            self._model_type = [
            #star
                    's---s-i', 's---smi', 
            #star + PLenvelope
                    's-p-smi', 's-p-hmi', 
            #star + PLenvelope + bipolar cavity
                    's-pbhmi', 's-pbsmi', 
            #star + disk 
                    'sp--h-i', 'sp--s-i',
                    'sp--hmi', 'sp--smi', 
            #star + ULenvelope
                    's-u-hmi', 's-u-smi',
            #star + ULenvelope + bipolar cavity
                    's-ubhmi', 's-ubsmi',  
            #star + disk + PLenvelope - Not present
            #
            #star + disk + ULenvelope
                    'spu-hmi', 'spu-smi', 
            #star + disk + ULenvelope + bipolar cavity
                    'spubhmi', 'spubsmi',  
                    ]
        else:
            self._model_type = model_type

        self._alpha = dict()
        self._wavelengths = [None,None]
        self._apertures   = [None]

    # how to deal with apertures?   
    def computeAlpha(self,color1,color2):
        for m in self._model_type:
            fname1=self._model_dir+m+'/convolved/'+color1+".fits"
            fname2=self._model_dir+m+'/convolved/'+color2+".fits"
            # convolved data cubes are 90000 fluxes by 20 apertures
            sed1 = ConvolvedFluxes.read(fname1)
            names = seds1.model_names
            self._apertures = seds1.apertures
            if np.sort(sed1.apertures) != np.sort(seds.apertures):
                raise Exception("apertures don't match")
            sed2 = ConvolvedFluxes.read(fname2)
            sed2.sort_to_match(seds1.model_names)
            self.wavelengths = [sed1.central_wavelegth.to("micron").value, 
                                sed2.central_wavelegth.to("micron").value]*u.micron

            # compute nu*Fnu for the two wavelengths
            sed2.nu = sed2.central_wavelength.to(u.Hz,equivalencies=u.spectral()) 
            sed1.nu = sed1.central_wavelength.to(u.Hz,equivalencies=u.spectral())
            self._alpha[m] =   np.log10(sed2.nu.value*sed2.flux.value) - np.log10(sed1.nu.value*sed1.flux.value)

           
class AlphaPlot():
    def __init__(self,wavelength,alpha,model=""):
        if len(alpha) != 2:
          raise Exception("Expected 2 alphas, got %d" % len(alpha))
        if len(wavelength) != 4:
          raise Exception("Expected 4 alphas, got %d" % len(wavelength))
        self._alpha = alpha
        self._wavelength = wavelength
        self._model = model

    def plot(self):
        ax, f = plt.figure()
        ax.scatter(self._alpha[0],self.alpha1[1],marker='.')
        ax[0].set_axislabel(r'$\alpha (%s - %s)'%(wavelength[1],wavelength[0]))
        ax[1].set_axislabel(r'$\alpha (%s - %s)'%(wavelength[3],wavelength[2]))
        plt.text("model: %s",self_model)
        plt.show()

if __name__ == "__main__":
       import filtermanage as fm
       a1 = AlphaMaker(model_type=['s---s-i'])
       w = []
       alph1 = a1.computeAlpha(fm.FORCAST_F371,fm.FORCAST_F242)
       w.extend(a1.wavelengths)
       alph2 = a1.computeAlpha(fm.IRAC4,fm.IRAC1)
       w.extend(a1.wavelengths)
       a = AlphaPlot(w,alpha1,alph2)
       a.plot()
    
