
import numpy as np

def radtran(wnm, od, tzl):
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


    nnu, nlyr = od.shape
    trans_lyr = np.vstack([np.ones((nnu)).T, np.exp(-od).T]).T
    transL = 1*trans_lyr
    for i in range(1, nlyr+1):
        transL[:,i] = transL[:,i-1]*trans_lyr[:,i]


    itop = len(tzl)-1
    if od.shape != (len(wnm), itop):
      raise ValueError('Input dimensions of optical depth are wrong.')
    
    # tave is not quite the average layer temperature.
    # Calculate it as the average temperature anyway
    tave = (tzl[:itop]+tzl[1:itop+1])/2
    
    # Downwelling radiance.  Start from surface and go up
    ilayers = list(range(itop))

    nulen  = len(wnm)
    rad    = np.zeros((nulen))

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


# .. A function to get the radiance for a model atmosphere (do not modify)
def get_my_radiance(co2, h2o, ch4, other, dT):
    
    nlyr_h2o = 9
    
    # .. Load in some files we will need
    nu = np.loadtxt('nu.txt')
    T0 = np.loadtxt('T.txt') - 10
    od_co2 = np.loadtxt('co2.txt')/367.71
    od_h2o = np.loadtxt('h2o.txt')/617.44
    od_ch4 = np.loadtxt('ch4.txt')/1.7
    od_other = np.loadtxt('other.txt')
    od_self = np.loadtxt('h2o_self.txt')/(617.44**2)
    
    # .. Get rid of some of the ringing in the od case
    # 781, 783, 786, 797, 800, 802
    inu = np.where(np.round(nu)==781)
    od_h2o[inu,:] = .000005
    inu = np.where(np.round(nu)==783)
    od_h2o[inu,:] = .000005
    inu = np.where(np.round(nu)==786)
    od_h2o[inu,:] = .000005
    inu = np.where(np.round(nu)==797)
    od_h2o[inu,:] = .000005
    inu = np.where(np.round(nu)==800)
    od_h2o[inu,:] = .000005
    inu = np.where(np.round(nu)==802)
    od_h2o[inu,:] = .000005
    inu = np.where(np.round(nu)==851)
    od_h2o[inu,:] = .000004

    
    # .. Create the new spectrum
    od = od_co2*co2 + od_ch4*ch4 + od_other*other
    od[:,:nlyr_h2o] = od[:,:nlyr_h2o] + od_h2o[:,:nlyr_h2o]*h2o \
    																	+ od_self[:,:nlyr_h2o]*h2o**2
    od[od<0] = 0
    T = T0 + dT
    
    rad = radtran(nu, od, T)
    
    my_legend = 'CO$_2$='+str(co2)+', H$_2$O='+str(h2o)+', CH$_4$='+str(ch4)+', other='+str(other)+'x'

    return nu, rad, my_legend


def plancknu(nu_icm,T):

  # function f = plancknu(nu_icm,T);
  #
  # spectral Planck function as function of wavenumbers (cm-1)
  #
  # [h]    = J*s
  # [c]    = m/s
  # [cbar] = cm/s
  # [k]    = J*K-1
  #
  #
  #    Note: LBLRTM uses ...
  #c    Constants from NIST 01/11/2002
  #
  #      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
  #     *     CLIGHT / 2.99792458E+10 /, 
  # 
 
  import numpy as np
   
  h    = 6.62606896e-34				# J s;  CODATA 2006
  c    = 2.99792458e8				# m/s;  NIST
  k    = 1.3806504e-23				# J K-1; CODATA 2006
  cbar = 100*c			  		# cm/s
  
  nu = nu_icm+0.
  
  # Ignore the possibility of input nu == 0 for now
  #if np.isscalar(nu):
  #  if nu==0:
  #    f = 0
  #else
  #  i0 = np.where(nu==0)[0]
  #  if len(i0) > 0:
  #    nu[i0] = 1e-6


  top    = 2 * h * cbar**3 * nu**3
  bottom =  c**2 *  ( np.exp(h*cbar*nu/(k*T))-1 )
  
  f      = cbar*top/bottom
 

  # now set f to zero for wavenumbers less than or equal to 0
  #if len(i0)>0:
  #  f[i0] = 0


  #[B]= cm*s-1*J*s*cm3*s-3*cm-3*m-2*s2
  #[B]= W m-2 cm

  return f

