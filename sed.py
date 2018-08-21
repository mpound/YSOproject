#!/usr/bin/env python
import numpy as np
from astropy import units as u
from astropy import coordinates as coord
from astropy.units.quantity import Quantity
from filtermanage import Photometry
import quantityhelpers as qh

class SED():
    """Simple spectral energy distribution class for input into SEDFitter"""
    def __init__(self,name,distance,disterr,ra=0*u.degree,dec=0*u.degree,av=None,averr=None):
       if not qh.isLength(distance):
          raise Exception("distance must be an Astropy Quantity with units length")
       if not qh.isQuantity(ra):
          raise Exception("Right Ascension must be an Astropy Quantity")
       if not qh.isQuantity(dec):
          raise Exception("Declination must be an Astropy Quantity")
       self._name     = str(name)
       self._distance = distance.to(u.pc)
       if type(disterr) == tuple:
           # error is +/-
           if disterr[1] > disterr[0]:
               self._disterr= (disterr[1].to(u.pc),disterr[0].to(u.pc))
           else:
               self._disterr= (disterr[0].to(u.pc),disterr[1].to(u.pc))
       else:
           self._disterr= (disterr.to(u.pc),-disterr.to(u.pc))
       if not qh.isLength(self._disterr[1]):
          raise Exception("minus distance error must be an Astropy Quantity with units length")
       if not qh.isLength(self._disterr[0]):
          raise Exception("plus distance error must be an Astropy Quantity with units length")
       if av >= 0:
           self._av     = av
       else:
           self._av     = None
       if averr >= 0:
           self._averr = averr
       else:
           self._averr = None
       self._coord = coord.ICRS(ra=ra,dec=dec)
       self._photometry = dict()

#    def fix_gaia_errors(self):
#        for b in ( "GAIA_G2" "GAIA_BP2" "GAIA_RP2" ):
#            if b in self._photometry:
#                if self._photometry[b].errormag()<0.004:
#                   self._photometry[b].errormag()<0.004:

    def set_upper_limits(self):
        for z in self._photometry.values():
            z.set_upper_limit(sn=3.0)

    def distrange(self):
        return [self._distance+self._disterr[1],self._distance+self._disterr[0]]

    def distrange_pc(self):
        return [self._distance.to(u.pc).value+self._disterr[1].to(u.pc).value,self._distance.to(u.pc).value+self._disterr[0].to(u.pc).value]

    def avrange(self):
        if self._av == None or self._averr == None: 
           return [None,None]
        if self._averr < self._av:
           return [self._av-self._averr,self._av+self._averr]
        else:
           return [0,self._av+self._averr]

    def addData(self,bandname,flux,error,validity,unit=None):
       if (np.ma.is_masked(flux) or np.ma.is_masked(error)) and validity != 0:
          print("flux and/or error is masked for Source %s band %s but validity is non-zero...skipping."%(self._name,bandname))
          return
       self._photometry[bandname.lower()] = Photometry(bandname,flux,error,validity,unit)

    def setvalidity(self,bandname,validity):
       """Reset the validity flag of an SED point"""
       self._photometry[bandname.lower()].setvalidity(validity)

    def header(self):
      """a descriptive header for human readability"""
      
      return "# d=%s                bands= %s " % (self._distance , self.bands())

    def bands(self):
       """Return the bands in this SED in the order sedfitterinput will return them"""
       line = " "
       for p in self._photometry.values():
           line += " %s" % p.band
       return line

#       w = ma.masked_array(np.zeros(len(self._photometry))
#       for i,s in enumerate(self._photometry.values(),0):
#           w[i] = s.wavelength.to(u.micron).value
#           if s.validity == 0:
#              w[i].mask == False

    def wavelengths(self):
       w = []
       for s in self._photometry.values():
           z = s.wavelength.to(u.micron).value
           w.append(z)
       return w*u.micron

    def fluxes(self):
       f = []
       for s in self._photometry.values():
           f.append(s.mjy())
       return f*u.mJy

    def sedfitterinput(self):
       """Return this SED as an input source for SEDFitter code"""
       #line = self.header() + "\n"
       # SEDfitter doesn't want spaces in names
       nospacename = self._name.replace(' ','_')
       line = "%s  %.3e  %.3e" % ( nospacename, self._coord.ra.degree, self._coord.dec.degree )
       photline = " "
       for p in self._photometry.values():
           photline+="%d "%p.validity
       for p in self._photometry.values():
           photline+="%5.4e %5.4e " % (p.mjy(), p.errormjy())
       
       return line+photline

if __name__ == "__main__":
       import filtermanage as fm
       pm = SED("MySource",100*u.pc,13.3*u.hour,-22*u.degree)
       # will raise exception
       try:
           pm = SED("testbaddistance",100*u.degree,13.3*u.hour,-22*u.degree)
       except Exception as e:
           print("EXCEPTION CAUGHT:",e)
       q = 1*u.Jy
       pm.addData(fm.SDSS_g,q/50, q/5000.0,0)
       pm.addData(fm.TWOMASS_K,q, q/100.0,1)
       # will raise warning
       pm.addData("3J",q/10, q/300.0,1)
       pm.addData(fm.SDSS_u,q/100, q/1000.0,0)
       print(pm.sedfitterinput())

       g = u.Magnitude(10)
       print(g)
       print(g.unit)
       print(qh.isMagnitude(g))
