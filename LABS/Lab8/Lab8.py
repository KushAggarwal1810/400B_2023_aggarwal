

# # Lab 8 : Star Formation 



import numpy as np
from astropy import units as u
from astropy import constants as const

import matplotlib
import matplotlib.pyplot as plt


# # Part A
# 
# Create a function that returns the SFR for a given luminosity (NUV, FUV, TIR, Halpha)
# 
# $Log( {\rm SFR} (M_\odot/year)) = Log(Lx (erg/s)) - Log(Cx)$ 
# 
# Including corrections for dust absorption 
# 
# Kennicutt & Evans 2012 ARA&A Equation 12 and Table 1, 2

def StarFormationRate(L,Type, TIR=0):
    """ thIS IS A FUNCTION THAT COMPUTES the star formation rate of a galaxy follwing the paper
    Kennicutt and evans 2012 eqn 12 (AR&A vol no:50)
    
    Inputs:
        
      L : Float
      luniosity of a g;aaxy in a given waveband which is ( erg/s)
      
      Type : 'string'
      
      the wavelength : 'FUV' ,'NUV', 'TIR', 'Halpha'
      
      TIR :float 
      
      total infrared luniosity in erg/s (defualt=0)
      
     OUTPUTS:
         
         SFR :float 
           Log of th star formation rate (Msun/yr)  
    """
    if (Type == 'FUV'):
        logCx = 43.35 # calibration from table 1 (K&E 2012)
        TIRc= 0.46  # correction for dust absorption from paper table 2(k and e 12)
    elif(Type == 'NUV'):
        logCx = 43.17
        TIRc= 0.27
    elif(Type == 'Halpha'):
        logCx = 41.27
        TIRc= 0.0024
        
    elif(Type == 'TIR'):
         logCx = 43.41
         TIRc= 0         # as we are already in infrared so no correction for infrared
    elif(Type == 'NUV'):
        print("Missing wavelength : FUV, NUV, Halpha, TIR")
    
    #Correct the luminoisty for dust using TIR
    Lnew = L + TIRc*TIR
    
    #star formation rate
    
    SFR= np.log10(Lnew)- logCx
    
    return SFR


# Let's try to reproduce SFRs derived for galaxies from UV luminosities measured with Galex. 
# 
# Using Table 1 from Lee et al. 2009
# https://ui.adsabs.harvard.edu/abs/2009ApJ...706..599L/abstract
# 
# We will use galaxy properties from NED:
# https://ned.ipac.caltech.edu/

# WLM Dwarf Irregular Galaxy
# From NED : WLM NUV Lumonisyty 1.71e7 Lsun
# From NED : WLM NIR luminosity 2.48e6 Lsun
# From NED : WLM FIR Luminosity 7.48e5 Lsun


NUV_WLM= 1.71e7

LsunErgS= const.L_sun.to(u.erg/u.s).value
print(LsunErgS)

NUV_WLM= 1.71e7*LsunErgS
TIR_WLM= 2.46e6*LsunErgS+ 7.84e5*LsunErgS

print(StarFormationRate(NUV_WLM, 'NUV',TIR_WLM))


#  WLM Dwarf Irregular Galaxy




#  N24 Sc galaxy


# # Part B Star formation main sequence
# 
# 1) Write a function that returns the average SFR of a galaxy at a given redshift. 
# 
# 2) What is the average SFR of a MW mass galaxy today? at z=1?
# 
# 3) Plot the SFR main sequence for a few different redshifts from 1e9 to 1e12 Msun.
# 
# 
# From Whitaker 2012:
# 
# log(SFR) = $\alpha(z)({\rm log}M_\ast - 10.5) + \beta(z)$
# 
# $\alpha(z) = 0.7 - 0.13z$
# 
# $\beta(z) = 0.38 + 1.14z - 0.19z^2$


# Step 1

def SFRMainSequence(Mstar,z):
    """ Function that computes the avergae sfr of a galaxy as a function of stellar mass
    
    INPUTS :
        
        mSTAR : 'FLOAT'
        STEALLR MASS OF THE GALXY IN SOLAR MASSES
        
        z: 'float'
        
        redshift
    OUTPUTS :
        
        logsfr : 'float'
        log(sfr(Msun/yr))
    """
    
    alpha= 0.7-0.13*z
    beta= 0.38 +1.14*z- 0.19*z**2
    logSFR= alpha*(np.log10(Mstar)-10.5)+beta
    
    return logSFR


    

# Step 2

# MW at z=0

MW_disk = 8e10
print(10**SFRMainSequence(MW_disk, 0))


# MW at z = 1

MW_disk = 8e10
print(10**SFRMainSequence(MW_disk, 1))  # this is avg sfr, to find total, integrate over the whole disk

# Step 3


# create an array of stellar masses

Mass= np.linspace(1e9,1e12)   # 112 massive elliptical galaxies, 1e9:sprial galxies





fig = plt.figure(figsize=(8,8), dpi=500)
ax = plt.subplot(111)

# add log log plots
plt.loglog(Mass,10**SFRMainSequence(Mass, 0), color= 'blue', linewidth=3, label= 'z=0')
plt.loglog(Mass,10**SFRMainSequence(Mass,5), color= 'red', linewidth=3, label= 'z=5')
plt.loglog(Mass,10**SFRMainSequence(Mass, 7.5), color= 'green', linewidth=3, label= 'z=7.5')
# plt.loglog(Mass,10**SFRMainSequence(Mass, 25), color= 'yellow', linewidth=3, label= 'z=25')
# plt.loglog(Mass,10**SFRMainSequence(Mass, 40), color= 'purple', linewidth=3, label= 'z=40')
plt.loglog(Mass,10**SFRMainSequence(Mass, 5), color= 'blue', linewidth=3, label= 'z=5')
plt.loglog(Mass,10**SFRMainSequence(Mass,8), color= 'red', linewidth=3, label= 'z=8')
plt.loglog(Mass,10**SFRMainSequence(Mass, 9), color= 'green', linewidth=3, label= 'z=9')
plt.loglog(Mass,10**SFRMainSequence(Mass, 6), color= 'yellow', linewidth=3, label= 'z=6')
plt.loglog(Mass,10**SFRMainSequence(Mass, 3), color= 'purple', linewidth=3, label= 'z=3')

# big z gives you nothing, as 

# Add axis labels
plt.xlabel('Log (Mstar (M$_\odot$))', fontsize=12)
plt.ylabel('Log(SFR (M$_\odot$/year))', fontsize=12)


#adjust tick label font size
label_size = 12
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

plt.show()

# # Part C  Starbursts

#Use your `StarFormationRate` code to determine the typical star formation rates for the following systems with the listed Total Infrared Luminosities (TIR): 

# Normal Galaxies: $10^{10}$ L$_\odot$

TIR_Normal=1e10*LsunErgS
print(10**StarFormationRate(TIR_Normal,'TIR'))
#
# LIRG: $10^{11}$ L$_\odot$
TIR_LIRG=1e11*LsunErgS
print(10**StarFormationRate(TIR_LIRG,'TIR'))
# ULIRG: $10^{12} $ L$_\odot$
# 
TIR_ULIRG=1e12*LsunErgS
print(10**StarFormationRate(TIR_ULIRG,'TIR'))
# HLIRG: $10^{13} $ L$_\odot$
TIR_HLIRG=1e13*LsunErgS
print(10**StarFormationRate(TIR_HLIRG,'TIR'))


# normal galaxies 




# LIRGs  




# ULIRGs




# HLIRGs






