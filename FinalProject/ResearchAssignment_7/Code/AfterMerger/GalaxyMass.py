# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 16:30:05 2023

@author: aggar
"""

import numpy as np
import astropy.units as u
from ReadFile import Read


def ComponentMass(filename,particletype):
    
    """
    This function computes the total mass associated to specific particle type for the input file
    
    Inputs :
        filename : 'txt file'
         The input file containing the data
         
         particle type : integer ( 1 = halo , 2 = buldge , 3 = buldge)
         
         This is the particle type for which the total mass is desired
         
    Outputs : 
        mass : 'astropy quantity'
        
        Total mass associated to the particle type input ( in 10^12 M_sun)
         
    """
    
    
    A= Read(filename)           # read function that we built in HW1.               
    data= A[2]         
    index = np.where(data['type']==particletype)[0]  
    mass =0
    for i in range(len(index)):
        a = data['m'][index[i]]
        mass = mass+a
    mass = np.round(mass,3)
    mass = (mass *1e10/1e12)*u.M_sun    # converting to 10^12 M_sun as given values are in 10^10 M_sun
    return mass