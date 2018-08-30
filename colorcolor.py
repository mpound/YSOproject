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
            names = sed1.model_names
            self._apertures = sed1.apertures
            sed2 = ConvolvedFluxes.read(fname2)
            sed2.sort_to_match(sed1.model_names)
            if not np.array_equal(sed1.apertures,sed2.apertures):
                raise Exception("apertures don't match")
            self.wavelengths = [sed1.central_wavelength.to("micron").value, 
                                sed2.central_wavelength.to("micron").value]*u.micron

            # compute nu*Fnu for the two wavelengths
            sed2.nu = sed2.central_wavelength.to(u.Hz,equivalencies=u.spectral()) 
            sed1.nu = sed1.central_wavelength.to(u.Hz,equivalencies=u.spectral())
            self._alpha[m] =   (np.log10(sed2.nu.value*sed2.flux.value) - np.log10(sed1.nu.value*sed1.flux.value))/(np.log10(sed2.nu.value)-np.log10(sed1.nu.value))

            return self._alpha[m]
           
class AlphaPlot():
    def __init__(self,wavelength,alpha,model=""):
        if len(alpha) != 2:
          raise Exception("Expected 2 alphas, got %d" % len(alpha))
        if len(wavelength) != 4:
          raise Exception("Expected 4 alphas, got %d" % len(wavelength))
        self._alpha = alpha
        print(np.shape(alpha))
        print(np.shape(np.transpose(self._alpha[0])[0]))
        self._wavelength = wavelength
        self._model = model

    def plot(self):
        a0 =np.transpose(self._alpha[0])[0]
        a1 =np.transpose(self._alpha[1])[0]
        plt.scatter(a0,a1,marker='.')
        plt.xlabel(r'$\alpha (%s - %s)'%(self._wavelength[1],self._wavelength[0]))
        plt.ylabel(r'$\alpha (%s - %s)'%(self._wavelength[3],self._wavelength[2]))
        plt.title("model: %s" % self._model)
        plt.show()

if __name__ == "__main__":
       import filtermanage as fm
       a1 = AlphaMaker(model_type=['s---s-i'])
       w = []
       alph1 = a1.computeAlpha(fm.FORCAST_F242,fm.FORCAST_F371)
       print(type(alph1))
       w.extend(a1.wavelengths)
       alph2 = a1.computeAlpha(fm.IRAC1,fm.IRAC4)
       print(type(alph2))
       print("A1 min,max,median,mean %.3f %.3f %.3f %.3f"%(np.min(alph1),np.max(alph1),np.median(alph1),np.mean(alph1)))
       print("A2 min,max,median,mean %.3f %.3f %.3f %.3f"%(np.min(alph2),np.max(alph2),np.median(alph2),np.mean(alph2)))
       w.extend(a1.wavelengths)
       a = AlphaPlot(w,[alph1,alph2],model="s---s-i")
       a.plot()
    
