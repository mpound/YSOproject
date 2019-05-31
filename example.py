#!/usr/bin/env python
import filtermanage as fm
import astropy.units as u

fsm = fm.FilterSetManager()
print("Filter sets imported:%s"%fsm.filtersetnames())

f = fsm.magtoflux(fm.SLOAN,fm.SDSS_u,10)
# note return value is Quantity
print(f)
print(f.to(u.Jy))
print(f.to(u.mJy))

m = fsm.fluxtomag(fm.SLOAN,fm.SDSS_u,156.85,mjy=True)
# magnitude Quantity 
print("Should be 10:",m)

# example using Quantity instead of scalar value input
q = 1000*u.mJy
m = fsm.fluxtomag(fm.SPITZER,fm.MIPS1,q)
print(m)
f = fsm.magtoflux(fm.SLOAN,fm.SDSS_u,0.0214)
print(f)

bandname = fm.SDSS_u
flux = 123.1*u.mJy 
error = 10*u.mJy
# See https://sedfitter.readthedocs.io/en/stable/data.html for explanation of
# validity.  
# validity 0 = invalid or unused, 1 = valid, 2 = lower limit, 3 = upper limit
validity = 3  


p = fm.Photometry(bandname,flux,error, validity)
print("flux is ",p.flux,"+/-",p.error,"[",p.error2,"]","[",p.error3,"]")
print("magnitude is ",p.magnitude,"+/-",p.errormag )
print("err val ",p.error.value)

bandname = fm.SDSS_g
flux = 10.1*u.mag
error = 1*u.mag
p = fm.Photometry(bandname,flux,error)
print("magnitude is ",p.magnitude,"+/-",p.errormag)
print("flux is ",p.flux,"+/-",p.error,"[",p.error2,"]","[",p.error3,"]")

bandname = fm.SDSS_g
flux = 10.1*u.mag
error = 0.1*u.mag
p = fm.Photometry(bandname,flux,error)
print("magnitude is ",p.magnitude,"+/-",p.errormag)
print("flux is ",p.flux,"+/-",p.error,"[",p.error2,"]","[",p.error3,"]")

error = 1.5*u.mag
p = fm.Photometry(bandname,flux,error)
print("magnitude is ",p.magnitude,"+/-",p.errormag)
print("flux is ",p.flux,"+/-",p.error,"[",p.error2,"]","[",p.error3,"]")

test = fm.Photometry(fm.SDSS_u,3*u.mag,1.*u.mag)
print("3 mag in sdss_u: ",test.flux)
test1 = fm.Photometry(fm.SDSS_u,test.flux,10.*u.mJy)
print("should be 3 mag: ",test1.magnitude)
test = fm.Photometry(fm.SDSS_u,9.194*u.mag,1.*u.mag)
print("9.194 mag in sdss_u: ",test.flux)
test1 = fm.Photometry(fm.SDSS_u,test.flux,10.*u.mJy)
print("should be 9.194 mag: ",test1.magnitude)

print("wavelength SDSS_u: ",test1.wavelength)
print("frequency SDSS_u: ",test1.frequency)

for error in [0.001, 0.01, 0.1, 1, 2 ]:
    p = fm.Photometry(fm.SDSS_u,10*u.mag,error*u.mag)
    print("{:>6.3f} {:>6.3f} {:>6.3f} {:>6.3f} {:>6.3f} {:>6.3f}".format( p.magnitude.value,p.errormag.value,p.flux.value,p.error.value,p.error2.value, p.error3.value))
    
