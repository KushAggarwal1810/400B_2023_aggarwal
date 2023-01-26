# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:55:43 2023

@author: aggar
"""

import numpy as np
import astropy.units as u
from ReadFile import Read     # importing the Read function that we made earlier


# this function prints the properties of any specific particle under study in the input file
   # at snaptime 0
# inputs : filename( in .txt format) , particletype(0 ,1 or 2) , particlenumber( integer)
# returns : mass of the particle ( in solar mass units) , 3d distance ( in light years)
            # 3d velocity ( in Km/s)
    
    
def ParticleInfo(filename,particletype,particlenumber):
    A= Read(filename)                       # reading the input file
    data= A[2]                              # Fetching the data array as that is the second element returned by Read
    index = np.where(data['type']==particletype)[0] # assigning index to particles with desired particle type : disc, halo etc.
    a = index[particlenumber-1]                # taing care of the fact that array start with 0.
    m = data['m'][a]                        # storing  mass 
    x = data['x'][a]                        # storing x distance
    y = data['y'][a]                        # storing y distance 
    z = data['z'][a]                        # storing z distance 
    vx = data['vx'][a]                      # storing  x velocity 
    vy = data['vy'][a]                      # storing x distance 
    vz = data['vz'][a]                      # storing x distance 
    m = m*1e10                              # converting mass to 1 solar mass 
    m= m*u.M_sun                            # assigning solar mass units to mass
    dist = np.sqrt(x**2+y**2 +z**2)*u.kpc    # calculating 3d distance and assigning kpc as units
    dist = dist.to(u.lyr)                  # converting kpc to light lyears
    velocity = np.sqrt(vx**2 +vy**2 +vz**2)*u.km /u. s      # calculating 3d velocity and assigning km/s as units.
    
    dist = np.around(dist,3)               # rounding of distance to 3 decimal places
   
    velocity = np.around(velocity,3)       # rounding of velocity to 3 decimal places 
    
    return m,dist,velocity                  # returns the value specified above in comments 

alpha = ParticleInfo('MW_000.txt', 2, 100)     # taking input from user about filename, particle type and particle number
print(f'the mass at snap number 0 is {alpha[0]:.3e}')
print(f'the 3D distance at snap number 0 is {alpha[1]:.3e}')


print(f'the 3D velocity at snap number 0 is {alpha[2]:.3e}')



    
    
    
     
    