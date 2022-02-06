# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 09:05:09 2015

@author: Penny Rowe and Steven Neshyba

Purpose:
    A simple illustration of inverse retrievals. Temperature and
    emissivity are retrieved for a grey-body cloud using optimal
    estimation. Atmospheric transmissivity is set to one.

Optimal estimation is described in detail by
    Rodgers, C. D., Inverse Methods for Atmospheric Sounding.
    World Scientific, 238 pp., 2000:

"""


import numpy as np
solve = np.linalg.solve 


# # # # # # # # # #     Plancknu     # # # # # # # # # # # # 
def plancknu(nu_icm_in,T):
  import copy  
  nu_icm = copy.deepcopy(nu_icm_in)
  #
  # spectral Planck function as function of wavenumbers (cm^-1)
  #
  # [h]    = J*s
  # [c]    = m/s
  # [cbar] = cm/s
  # [k]    = J*K-1
  # [B]    = cm*s-1*J*s*cm3*s-3*cm-3*m-2*s2
  # [B]    = W m-2 cm
  #
  h    = 6.62606896e-34				# J s;  CODATA 2006
  c    = 2.99792458e8				# m/s;  NIST
  k    = 1.3806504e-23				# J K-1; CODATA 2006
  cbar = 100*c			  		    # cm/s

  indzero = np.where(nu_icm==0)            # avoid divide-by-zero  
  nu_icm[indzero]=.1
  
  top        = 2 * h * cbar**3 * nu_icm**3
  bottom     =  c**2*  ( np.exp(h*cbar*nu_icm/(k*T))-1 )
  f          = cbar*top/bottom 
  f[indzero] = 0


  return f
  
# # # # # # # # # #     The forward model     # # # # # # # # # # # # 
def greybody(nu,X):
    T   = X[0]
    Eps = X[1]
    
    y = 1e3*plancknu(nu,T)  * Eps
    
    return y
  
  
  
# # # # # # # # # #     Main Code     # # # # # # # # # # # # 


# .. Set true cloud temperature and emissivity
#    1 RU = 1 mW / (m2 sr cm-1)
T = 300.                  # Temperature, K
Eps = .5
ErrLevelEst    = 0.1      # Estimate of noise level, RU
NoiseLevelTrue = 0.1      # Noise (random),RU
Nobs = 5                  # Number of observations, or wavenumbers
bias = 0                  # Bias error

# .. wavenumbers
nu = np.linspace(500, 1500, Nobs)

# .. A priori and statistics: 
#   note: X is temperature, emissivity, y is radiance at nu
Xa     = np.array([273., .8])              # a priori X
Xfg    = np.array([273., .8])              # first guess X
Sa_vec = np.array([30**2., 1])             # variance for Xa
Se_vec = ErrLevelEst**2 * np.ones((Nobs))  # variance in measurement, y
Niters = 5                                 # number of iterations

# .. Make the observed spectrum, y
noise = NoiseLevelTrue*np.random.randn(Nobs)
yobs  = (1e3*plancknu(nu,T) + noise + bias) * Eps





# .. Set variables that won't change
Sa     = np.diag(Sa_vec)
Se     = np.diag(Se_vec)
inv_Sa = np.diag(1/Sa_vec);
yfg    = greybody(nu,Xfg)
yn_1   = yfg + 0.
Xn_1   = Xfg + 0.
dT     = .1
dEps   = .01
Kn     = np.ones((Nobs,2))  
 


# .. Do the retrieval

for iter in range(Niters):

    # .. Get kernels. Xp is the perturbed Xn_1
    #    Temperature, T
    Xp      = Xn_1 + 0.; Xp[0] += dT
    yp      = greybody(nu,Xp)
    Kn[:,0] = (yp - yn_1 ) / dT     # Temperature

    #     Emissivity, Eps
    Xp = Xn_1 + 0; Xp[1] += dEps
    yp = greybody(nu,Xp)
    Kn[:,1] = (yp - yn_1 ) / dEps
    
    
    
    # .. Invert as in Rodgers eqn. 5.9 (n-form)
    KT_Sem1   = ( solve(Se.T, Kn) ).T
    KT_Sem1_K = np.dot( KT_Sem1 , Kn )
      
    term2 = inv_Sa + KT_Sem1_K 
    term3 = np.dot( KT_Sem1 , ( yobs-yn_1 + np.dot(Kn,(Xn_1-Xa))  ) )
    Xn    = Xa + solve(term2, term3)
    
    yn    = greybody(nu,Xn)
    
    # .. Set up for next iteration
    Xn_1 = Xn + 0.
    yn_1 = yn + 0.

    print(Xn)