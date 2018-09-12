#!/usr/bin/env python 
# read clasxx.data and write to ipac format
#File Summary:
#--------------------------------------------------------------------------------
# FileName      Lrecl  Records   Explanations
#--------------------------------------------------------------------------------
#ReadMe            80        .   This file
#clasi-ii.dat     354   133980   Class I/II YSO candidates
#clasiii.dat      354   608606   Class III and more evolved YSO candidates
#--------------------------------------------------------------------------------
#
#Byte-by-byte Description of file: clasi-ii.dat clasiii.dat
#--------------------------------------------------------------------------------
#   Bytes Format Units   Label     Explanations
#--------------------------------------------------------------------------------
#   1- 20  A20   ---     Name      AllWISE Identifier
#  22- 32  F11.7 deg     RAdeg     Right ascension (J2000.0)
#  34- 44  F11.7 deg     DEdeg     Declination (J2000.0)
#  46- 56  F11.7 deg     GLON      Galactic longitude
#  58- 68  F11.7 deg     GLAT      Galactic latitude
#  70- 75  F6.3  mag     W1mag     Instrumental profile-fit photometry
#                                   magnitude, band 1
#  77- 81  F5.3  mag   e_W1mag     Instrumental profile-fit photometry flux
#                                   uncertainty in mag units, band 1
#  83- 87  F5.1  ---     W1snr     Instrumental profile-fit photometry
#                                   S/N ratio, band 1
#  89-100  F12.7 ---     w1rchi2   Instrumental profile-fit photometry reduced
#                                   chi^2^, band 1
# 102-107  F6.3  mag     W2mag     Instrumental profile-fit photometry
#                                   magnitude, band
# 109-113  F5.3  mag   e_W2mag     Instrumental profile-fit photometry
#                                   flux uncertainty in mag units, band 2
# 115-120  F6.1  ---     W2snr     Instrumental profile-fit photometry
#                                   S/N ratio, band 2
# 122-133  F12.7 ---     W2rchi2   Instrumental profile-fit photometry
#                                   reduced chi^2, band 2
# 135-140  F6.3  mag     W3mag     Instrumental profile-fit photometry
#                                   magnitude, band 3
# 142-146  F5.3  mag   e_W3mag     Instrumental profile-fit photometry
#                                   flux uncertainty in mag units, band 3
# 148-152  F5.1  ---     W3snr     Instrumental profile-fit photometry
#                                   S/N ratio, band 3
# 154-165  F12.7 ---     W3rchi2   Instrumental profile-fit photometry
#                                   reduced chi^2, band 3
# 167-172  F6.3  mag     W4mag     Instrumental profile-fit photometry
#                                   magnitude, band 4
# 174-178  F5.3  mag   e_W4mag     Instrumental profile-fit photometry
#                                   flux uncertainty in mag units, band 4
# 180-185  F6.1  ---     W4snr     Instrumental profile-fit photometry
#                                   S/N ratio, band 4
# 187-198  F12.7 ---     W4rchi2   ?=- Instrumental profile-fit photometry
#                                   reduced chi^2, band 4
# 200-211  F12.7 ---     rchi2     ?=- Instrumental profile-fit photometry
#                                   reduced chi squared, total
#     213  I1    ---     extflag   Probability that source morphology is not
#                                       consistent with single PSF
# 215-218  A4    ---     phQual    Photometric quality of each band
#                                   (A=highest, U=upper limit)
# 220-222  I3    ---     W1nm      Number of profile-fit flux measurements for
#                                   source with SNR >= 3, band 1
# 224-226  I3    ---     W1m       Number of profile-fit flux measurements for
#                                   source, band 1
# 228-230  I3    ---     W2nm      Number of profile-fit flux measurements for
#                                   source with SNR >= 3, band 2
# 232-234  I3    ---     W2m       Number of profile-fit flux measurements for
#                                   source, band 2
# 236-238  I3    ---     W3nm      Number of profile-fit flux measurements for
#                                   source with SNR >= 3, band 3
# 240-242  I3    ---     W3m       Number of profile-fit flux measurements for
#                                   source, band 3
# 244-246  I3    ---     W4nm      Number of profile-fit flux measurements for
#                                   source with SNR >= 3, band
# 248-250  I3    ---     W4m       Number of profile-fit flux measurements for
#                                   source, band 4
# 252-261  I10   ---     2MASSKey  Closest associated 2MASS All-Sky Release
#                                   PSC key
# 263-268  F6.3  mag     Jmag      J magnitude entry of the associated
#                                   2MASS All-Sky PSC source
# 270-274  F5.3  mag   e_Jmag      J photometric uncertainty of the associated
#                                   2MASS All-Sky PSC source
# 276-281  F6.3  mag     Hmag      H magnitude entry of the associated
#                                   2MASS All-Sky PSC source
# 283-287  F5.3  mag   e_Hmag      H photometric uncertainty of the associated
#                                   2MASS All-Sky PSC source
# 289-294  F6.3  mag     Ksmag     Ks magnitude entry of the associated
#                                   2MASS All-Sky PSC source
# 296-300  F5.3  mag   e_Ksmag     Ks photometric uncertainty of the associated
#                                   2MASS All-Sky PSC source
# 302-334  A33   ---     SName     Source identifier as listed in SIMBAD (1)
# 336-338  A3    ---     SType     SIMBAD maintype of the identified source (1)
# 340-344  F5.3  arcsec  SDist     ?=0 Angular distance between AllWISE position
#                                   and the SIMBAD counterpart (1)
# 346-354  F9.7  ---     tau       Optical depth at 353GHz from Planck
#                                   HFI_CompMap_ThermalDustModel_2048_R1.20.fits
#--------------------------------------------------------------------------------

slices = [
  [1, 20],
  [22, 32],
  [34, 44],
  [46, 56],
  [58, 68],
  [70, 75],
  [77, 81],
  [83, 87],
  [89,100],
 [102,107],
 [109,113],
 [115,120],
 [122,133],
 [135,140],
 [142,146],
 [148,152],
 [154,165],
 [167,172],
 [174,178],
 [180,185],
 [187,198],
 [200,211],
 [213,214],
 [215,218],
 [220,222],
 [224,226],
 [228,230],
 [232,234],
 [236,238],
 [240,242],
 [244,246],
 [248,250],
 [252,261],
 [263,268],
 [270,274],
 [276,281],
 [283,287],
 [289,294],
 [296,300],
 [302,334],
 [336,338],
 [340,344],
 [346,354]
]

# need to be zero based

for j in slices:
    j[0] = j[0]-1
    j[1] = j[1]-1
#key = variable name, value = unit
cols = {
   "Name":"",
   "RAdeg":"degree",   
   "DEdeg":"degree",   
   "GLON":"degree",    
   "GLAT":"degree",    
   "W1mag":"mag",   
   "e_W1mag":"mag",   
   "W1snr":"",   
   "w1rchi2":"", 
   "W2mag":"mag",   
   "e_W2mag":"mag",   
   "W2snr":"",   
   "W2rchi2":"", 
   "W3mag":"mag",   
  "e_W3mag":"mag",   
   "W3snr":"",   
   "W3rchi2":"", 
   "W4mag":"mag",   
  "e_W4mag":"mag",   
   "W4snr":"",   
   "W4rchi2":"", 
   "rchi2":"",   
   "extflag":"", 
   "phQual":"",  
   "W1nm":"",    
   "W1m":"",     
   "W2nm":"",    
   "W2m":"",     
   "W3nm":"",    
   "W3m":"",     
   "W4nm":"",    
   "W4m":"",     
   "2MASSKey":"",
   "Jmag":"mag",    
   "e_Jmag":"mag",    
   "Hmag":"mag",    
   "e_Hmag":"mag",    
   "Ksmag":"mag",   
   "e_Ksmag":"mag",   
   "SName":"",   
   "SType":"",   
   "SDist":"arcsec",   
   "tau":""
}

i = 0 
head_format = []
var_format = []
print("|",end="")
for c in cols:
    # sometimes variable name or unit lengths are larger than the slice
    num = max(slices[i][1]-slices[i][0],len(c))
    num = max(num,len(cols[c]))
    head_format.append('{:>%d} | '%(num))
    var_format.append('{:>%d}   '%(num))
    print(head_format[i].format(c),end="")
    i=i+1
# NOT print("\n") because that is TWO newlines. Unless
# we do print("\n",end="").  gotta love python
print("")

# data type
i = 0 
print("|",end="")
for c in cols.values():
    if c != "":
        print(head_format[i].format("float"),end="")
    else:
        print(head_format[i].format("char"),end="")
    i=i+1

# data unit
i = 0 
print("|",end="")
for c in cols.values():
    print(head_format[i].format(c),end="")
    i=i+1

print("")
with open("clasi-ii.dat") as f:
#with open("clasiii.dat") as f:
    for line in f:
        vals = [line[slice(*slc)] for slc in slices]
        i=0
        print(" ",end="")
        for v in vals:
            #if v=="?": v="null"
            print(var_format[i].format(v),end="")
            i=i+1
        print("")

print("")

