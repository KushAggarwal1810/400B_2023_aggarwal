
# # In Class Lab 1
# Must be uploaded to your Github repository under a "Labs/Lab1" folder by 5 PM on Jan 31st 2023

# ## Part A:  The Local Standard of Rest
# Proper motion of Sgr A* from Reid & Brunthaler 2004
# $\mu = 6.379$ mas/yr 
# 
# Peculiar motion of the sun, $v_\odot$ = 12.24 km/s  (Schonrich 2010)
# 
# 
# $v_{tan} = 4.74 \frac{\mu}{\rm mas/yr} \frac{R_o}{\rm kpc} = V_{LSR} + v_\odot$
# 
# 
# ### a)
# 
# Create a function called VLSR to compute the local standard of res (V$_{LSR}$).
# 
# The function should take as input: the solar radius (R$_o$), the proper motion (mu)
# and the peculiar motion of the sun in the $v_\odot$ direction.
# 
# Compute V$_{LSR}$ using three different values R$_o$: 
# 1. Water Maser Distance for the Sun :  R$_o$ = 8.34 kpc   (Reid 2014 ApJ 783) 
# 2. GRAVITY Collaboration Distance for the Sun:  R$_o$ = 8.178 kpc   (Abuter+2019 A&A 625)
# 3. Value for Distance to Sun listed in Sparke & Gallagher : R$_o$ = 7.9 kpc 
# 



# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from astropy import constants as const # import astropy constants



def VLSR(Ro, mu= 6.379,vsun= 12.24*u.km/u.s): # fixing mu and vsun
   """ This function will compute the velocity at local standard of rest.
   VLSR=  4.74*mu*Ro - vsun
   Inputs :
       Ro : 'astropy quantity'
            The distance from the sun to the galactic center in kpc.
            
       mu  : 'float'
            The proper motion of Sag A* in mas/yr .
            Default is from Reid & Brunthaler 2004
    vsun   : 'astropy quantity'
             The peculiar motion of Sun in v direction ( km/s) 
             Default is from Schonrich + 2010
             
    Outputs:
         VLSR : 'Astropy Quantity'
          The velocity if local standard of rest ( km/s)
   """
   
   return 4.74*mu*(Ro/u.kpc)*u.km/u.s - vsun
    
    
# define our distances tot he galactic center from the sun

RoReid = 8.34 *u.kpc # Distance from Reid et al . 2014 in kpc
RoGravity = 8.178*u.kpc # distance from the Gravity Collab Abuter + 2019 in kpc
RoSG = 7.9 *u.kpc # Distance from the textbook Sparke and Gallaghar


# Compute VLSR using Ro  from Reid 2014

VLSR_Reid = VLSR(RoReid)
print(VLSR_Reid)

# Compute VLSR using Ro from gravity collab

VLSR_Gravity = VLSR(RoGravity)
print(np.round(VLSR_Gravity))

#Compute VLSR using Ro from Sparke and Gallaghar

VLSR_SG = VLSR(RoSG)
print(np.round(VLSR_SG))



# ### b)
# 
# compute the orbital period of the sun in Gyr using R$_o$ from the GRAVITY Collaboration (assume circular orbit)
# 
# Note that 1 km/s $\sim$ 1kpc/Gyr


def TorbSun(R,V):
    
    """ This function will compute the orbital period of the sun
     T= 2 pi R / V
     
     Inputs : R 'Astropy Quantity'
                Distance in kpc ( Distance to the Galactic center)
                
             V : 'Asropy Quantity'
                    Velocity in Km/s ( Velocity of the sun in v direction)
    Outputs : 
         'astropy quantity'
         Orbital Period in GYr
                    
    """ 
    VkpcGyr = V. to(u.kpc/u.Gyr) # Converting v from km/s to kpc/Gyr
    T = 2*np.pi*R/VkpcGyr
    return T

# Velocity of the sun  = VLSR +peculair motion
VsunPeculiar = 12.24*u.km/u.s
Vsun = VLSR_Gravity +VsunPeculiar
#Compute the orbital period of sun
# use R0 from GRavity collaboration  
  
T_Grav = TorbSun(RoGravity,Vsun)
print(T_Grav)




# ### c)
# 
# Compute the number of rotations about the GC over the age of the universe (13.8 Gyr)

# Age of the Universe/Orbital Period

Age = 13.8 *u.Gyr # Age of the universe
print(Age/T_Grav)




# ## Part B  Dark Matter Density Profiles
# 
# ### a)
# Try out Fitting Rotation Curves 
# [here](http://wittman.physics.ucdavis.edu/Animations/RotationCurve/GalacticRotation.html)

""" 0, 86, 11"""
# 
# ### b)
# 
# 
# In the Isothermal Sphere model, what is the mass enclosed within the solar radius (R$_o$) in units of M$_\odot$? 
# 
# Recall that for the Isothermal sphere :
# $\rho(r) = \frac{V_{LSR}^2}{4\pi G r^2}$
# 
# Where $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$, r is in kpc and $V_{LSR}$ is in km/s
# 
# What about at 260 kpc (in units of  M$_\odot$) ? 

# density profile rho = VLSR^2/(4*pi*G*R^2)
# Mass = Integrate rho dV
#  = rho*4*pi*r^2 dr
# = VLSR^2/(4*pi*G*r^2)*4*pi*r^2 dr
# VLSR ^2/ G*r 
# Gravitational constant 

Grav = const. G.to(u.kpc**3/u.Gyr**2/u.Msun)
def MassIso(r, VLSR):
    """ This function will compute the dark matter mass enclosed in a given distnace assuming a isotherma
    sphere model for th drak matter
    M = VLSR**2/  G*r
    Inputs :
        r : "astropy units"
           distance to the galactic center in kpc
        VLSR: "astropy qunatity"
        
          velcoity of locla standard of rest in km/s
    Outputs :
        M : Mass encosed within r in units of Msun
        
    """
   
    VLSRkpcGyr = VLSR.to(u.kpc/u.Gyr) # converting km/s to kpc/Gyr
    M = VLSRkpcGyr**2 / Grav*r     # Mass for isothermal sphere
    return M

MIsoSolar= MassIso(RoGravity, VLSR_Gravity)
print(MIsoSolar)

# Print in Scientific Notation

print(f"{MIsoSolar:.2e}")

# compute mass within 260kpc
MIso260 =  MassIso(260*u.kpc , VLSR_Gravity)
print(f"{MIso260:.2e}")


# ## c) 
#
# The Leo I satellite is one of the fastest moving satellite galaxies we know. 
# 
# 
# It is moving with 3D velocity of magnitude: Vtot = 196 km/s at a distance of 260 kpc (Sohn 2013 ApJ 768)
# 
# If we assume that Leo I is moving at the escape speed:
# 
# $v_{esc}^2 = 2|\Phi| = 2 \int G \frac{\rho(r)}{r}dV $ 
# 
# and assuming the Milky Way is well modeled by a Hernquist Sphere with a scale radius of $a$= 30 kpc, what is the minimum mass of the Milky Way (in units of M$_\odot$) ?  
# 
# How does this compare to estimates of the mass assuming the Isothermal Sphere model at 260 kpc (from your answer above)


# Potential for a hernquist profile 
# Phi = -G*M/(r+a)

# using the potential for a hernquist pofle the equation for thr escape seed becomes:
# vesc**2 = 2*G*M/ (r+a)

# Reaarange the escape speed equation for M , what I want is mass

# So , M = vesc**2 / 2/ G *(r+a)


def MassFromVesc(vesc ,r,a):
    """ this function determines hte total mass needed for a given escape speed 
    assuimng a hernquist prifle for the dark matter halo
    M = vesc**2*(r+a) /2/G
    
    Inputs : 
        vesc : 'astropy quantity'
        The escape speed in km/s ( or speed of the satellite)
        r : 'astropy qty'
          The distance from galactic center ( in kpc)
        
        a : 'astropy qty'
        The Hernquist scale length (kpc)
        
    Outputs : 
        M : 'astropy quantity'
        The total mass within r in Msun
        
    """
    
    vescKpcGyr= vesc.to(u.kpc/u.Gyr) # converting velcoity units to kpc/Gyr
    M= vescKpcGyr**2/2/Grav*(r+a) # required mass
    return M

VLeoI = 196*u.km/u.s # Speed of Leo 1
a = 30 *u.kpc # Scale Radius for the Hernquist Halo
r= 260*u.kpc #  Galactic center distance of Leo 1

# Compute the mass needed to keep Leo1 bound

MLeoI = MassFromVesc(VLeoI, r, a)

print(f"{MLeoI:.2e}")

print(MIso260/MLeoI)  # remember mass increase linearly with radius





