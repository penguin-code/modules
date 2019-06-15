
import numpy as np
from plancknu import plancknu




#function [rad,Rbelow,Rabovecloud]=radtran_lorestrans4fast_ud(...
#  wnm,od,tzl,tave,izeroout,bltype,cloudlayer,direction,radsurf)

def radtran_t(wnm, transL, tzl):
#
# function [rad,Rbelow,Rabovecloud]=radtran_lorestrans4fast_ud(...
#  wnm,od,tzl,tave,izeroout,bltype,cloudlayer,direction,radsurf)
#
# by Penny Rowe
# May 30, 2001
# Updated Aug. 7, 2002
# Updated Feb. 2012: added odDir
# Updated May 23, 2013: 
#   Improved comments
#   Changed rdODsingle to rdODprec, added "prec" (precision) to inputs
#   Changed name from radtranC to radtran_lorestrans
#   Added possibilities "lowClough13" and "lowCloughPade" to bltype
#   No longer cd to directory where ods are
# Updated Jan. 11, 2013
#   Can do up or downwnwelling radiance, based on direction
#
#
# Purpose: 
#   Get low- or perfect-resolution radiance (and corresponding wavenumber)
#   by performing radiative transfer in the same way as LBLRTM,
#   except that Eqn (13) of Clough et al. 1992 is used instead of
#   getting the Planck function empirically using a Pade coefficient
#   as in Eqn (15), *or* using the Pade function, or using the method
#   of Wiscombe 1976.  Also, can do low or perfect resolution.  These
#   Choice depend on input variable bltype.
#
#   Low-resolution is acheived by reducing the resolution of the
#   transmittances, not the radiances
#
# .. Inputs to function call:
#    wnm  = wavenumber, numpy 1-d array
#    od   = layer optical depths, numpy array [itop x len(wnm)]
#    itop = number of layers for which there are ods
#    tzl  = temperature at layer bottom, from surface up
#    tave = "ave" temperature, as defined by LBLRTM
#    viewangle = viewing angle, in degrees
#
#    Old, obsolete:
#    dnu  = desired resolution
#    bltype: lowClough13, lowClough/lowCloughPade, hiClough13, lowWiscombe
#
#    izeroout  removed!
#
# .. Outputs:
#    rad: radiance, numpy 1-d array
#
# .. Updates:
#    Accounting for off-zenith viewing angle. Done by dividing
#    the optical depth by cosine(angle_in_radians) throughout.
# # # # # # # # # # # #  # # # # #  # # # #  # # # #  # # # # # # #

		print(transL.shape)
    nnu, nlyr = transL.shape
    #trans_lyr = np.vstack([np.ones((nnu)).T, np.exp(-od).T]).T
    #transL = 1*trans_lyr
    # Watch out, I think t(0)=1, so we need to omit the first layer!
    od = 1*transL
    od[:,0] = -log(transL[:,0])
    for i in range(1, nlyr+1):
        od[:,i] = -log(transL[:,i]) - od[:,i-1] ;
    #    transL[:,i] = transL[:,i-1]*trans_lyr[:,i]


    itop = len(tzl)-1
    if transL.shape != (len(wnm),itop):
      raise ValueError('Input dimensions of trans are wrong.')
    
    # tave is not quite the average layer temperature.
    # Calculate it as the average temperature anyway
    tave = (tzl[:itop]+tzl[1:itop+1])/2
    
    # Downwelling radiance.  Start from surface and go up
    ilayers = list(range(itop))

    nulen = len(wnm)
    rad = np.zeros((nulen))

    # get transmittance, rad., from the bottom layer up
    for ilayer in ilayers:
        
      # Uses low resolution layer optical depth
      # See Clough et al. 1992 Eqn. (13)
      # PMR May 23, 2013
      Bb   =  plancknu(wnm, tave[ilayer])
      Bu   =  plancknu(wnm, tzl[ilayer])
      
      atau = 0.278*od[:,ilayer]                       # Pade coeff. x tau
      BL   = (Bb + atau*Bu)* (1+atau)**-1 ;                    # Eqn (15)
      rad0 = -BL * (transL[:,ilayer+1]-transL[:,ilayer]) ;     # Eqn (14)
      rad += rad0

    # QC
    if sum(np.imag(rad))<1e-10:
      rad=np.real(rad) 
    else:
      raise ValueError('Non-trivial imaginary part')
    
    rad     = 1e3*rad

    return rad
