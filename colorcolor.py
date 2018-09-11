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
        self._nufnu  = dict()

    # how to deal with apertures?   
    def computeNuFnu(self,color1):
        for m in self._model_type:
            fname1=self._model_dir+m+'/convolved/'+color1+".fits"
            # convolved data cubes are 90000 fluxes by 20 apertures
            sed1 = ConvolvedFluxes.read(fname1)
            names = sed1.model_names
            self._apertures = sed1.apertures
            # compute nu*Fnu for the wavelengths
            sed1.nu = sed1.central_wavelength.to(u.Hz,equivalencies=u.spectral())
            self._nufnu[m] =  np.log10(sed1.nu.value*sed1.flux.to("Jy")*1.0E26) 
            print("fnu ",len(self._nufnu))
        return self._nufnu

    def computeAlpha(self,band1,band2):
        """Compute dLog(lambda*Flambda)/dLog(lambda) where nu is the frequency and Fnu is the flux 
           density.  Specifically, this function will return
               (log10(lambda1*flambda1)-log10(lambda22*flambda2))/(log10(lambda1)-log10(lambda2))
           where lambda1/2,flambda1/2 correspond to input band1/2 respectively.

           Parameters:
               band1 - string indicating input band 1. Must be defined in filtermanage
               band2 - string indicating input band 2. Must be defined in filtermanage
        """
        for m in self._model_type:
            fname1=self._model_dir+m+'/convolved/'+band1+".fits"
            fname2=self._model_dir+m+'/convolved/'+band2+".fits"
            # convolved data cubes are 90000 fluxes by 20 apertures
            sed1 = ConvolvedFluxes.read(fname1)
            names = sed1.model_names
            self._apertures = sed1.apertures
            sed2 = ConvolvedFluxes.read(fname2)
            sed2.sort_to_match(sed1.model_names)
            if not np.array_equal(sed1.apertures,sed2.apertures):
                raise Exception("apertures don't match")
            self.wavelengths = [sed1.central_wavelength.to("meter").value, 
                                sed2.central_wavelength.to("meter").value]*u.meter

            # compute nu*Fnu for the two wavelengths
            sed1.nu = sed1.central_wavelength.to(u.Hz,equivalencies=u.spectral())
            sed2.nu = sed2.central_wavelength.to(u.Hz,equivalencies=u.spectral()) 
            # flux is in mJy
            print("flux %s"%sed1.flux[0].unit.decompose())
            sed1.nufnu = sed1.central_wavelength.to("m")*sed1.flux.to("Jy")*1.0E26
            sed2.nufnu = sed2.central_wavelength.to("m")*sed2.flux.to("Jy")*1.0E26
            print("flux %s"%sed1.nufnu.unit.decompose())
            diff = np.log10(sed1.nufnu.value)-np.log10(sed2.nufnu.value)
            print(diff)
            #self._alpha[m] =   diff/(sed1.nu.value-sed2.nu.value)
            self._alpha[m] =   diff/(np.log10(self.wavelengths[0].value)-np.log10(self.wavelengths[1].value))
            #print(self._alpha[m])

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
        a0 =np.transpose(self._alpha[0])[-1]
        a1 =np.transpose(self._alpha[1])[-1]
        plt.scatter(a0,a1,marker='.')
        plt.xlabel(r'$\alpha$ (%s - %s)'%(self._wavelength[0].to("micron"),self._wavelength[1].to("micron")))
        plt.ylabel(r'$\alpha$ (%s - %s)'%(self._wavelength[2].to("micron"),self._wavelength[3].to("micron")))
        plt.title("model: %s" % self._model)
        plt.show()

if __name__ == "__main__":
       import filtermanage as fm
       mymodel = 's---s-i'
       #mymodel = 'sp--h-i'
       #mymodel = 's-p-smi'
       #mymodel = 'spubhmi'
       a1 = AlphaMaker(model_type=[mymodel])
       w = []
       alph1 = a1.computeAlpha(fm.FORCAST_F113,fm.FORCAST_F371)
       #print(type(alph1))
       w.extend(a1.wavelengths)
       alph2 = a1.computeAlpha(fm.IRAC1,fm.IRAC4)
       #nf1   = a1.computeNuFnu(fm.IRAC1)
       #nf2   = a1.computeNuFnu(fm.FORCAST_F371)
       print(type(alph2))
       print("A1 min,max,median,mean %.3f %.3f %.3f %.3f"%(np.min(alph1),np.max(alph1),np.median(alph1),np.mean(alph1)))
       print("A2 min,max,median,mean %.3f %.3f %.3f %.3f"%(np.min(alph2),np.max(alph2),np.median(alph2),np.mean(alph2)))
       w.extend(a1.wavelengths)
       if True:
           a = AlphaPlot(w,[alph1,alph2],model=mymodel)
       else:
           a = AlphaPlot(w,[nf1,nf2],model=mymodel)
       a.plot()
    
