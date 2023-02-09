# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 05:39:49 2023

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

"""
Priniting the total type 1 mass in each galaxy
"""
A = ComponentMass('M31_000.txt',1)
A= np.round(A,3)
print("Mass of type 1 particles in M31 is", A)

B = ComponentMass('M33_000.txt',1)
B= np.round(B,3)
print("Mass of type 1 particles in M33 is", B)

C= ComponentMass('MW_000.txt',1)
C= np.round(C,3)
print("Mass of type 1 particles in MW is", C)

"""
Priniting the total type 2 mass in each galaxy
"""
D = ComponentMass('M31_000.txt',2)
D= np.round(D,3)
print("Mass of type 2 particles in M31 is", D)

E = ComponentMass('M33_000.txt',2)
E= np.round(E,3)
print("Mass of type 2 particles in M33 is", E)


F= ComponentMass('MW_000.txt',2)
F= np.round(F,3)

print("Mass of type 2 particles in MW is", F)

"""
Priniting the total type 3 mass in MW and M31
"""
G= ComponentMass('M31_000.txt',3)
G= np.round(G,3)
print("Mass of type 3 particles in M31 is", G)

H= ComponentMass('MW_000.txt',3)
H= np.round(H,3)
print("Mass of type 3 particles in MW is", H)

print("Mass of type 3 particles in M33 is", 0)

"""
Printing total mass of each galaxy
"""
MW_totalmass = C+F+H
MW_totalmass = np.round(MW_totalmass,3)
print("Total mass of MW is", MW_totalmass)
M31_Totalmass= A+D+G
M31_Totalmass = np.round(M31_Totalmass,3)
print("Total mass of M31 is", M31_Totalmass)

M33_Totalmass = B+E
M33_Totalmass= np.round(M33_Totalmass,3)

print("Total mass of M33 is",M33_Totalmass)

M_Localgroup = MW_totalmass +M31_Totalmass+M33_Totalmass

M_Localgroup = np.round(M_Localgroup,3)

print("Total mass of the local group is",M_Localgroup)

"""
Calculating fbar : the ratio of stellar mass to total mass ( dark +matter) for each gaalaxy
"""
fbar_MW= (F+H)/MW_totalmass
fbar_MW = np.round(fbar_MW,3)
print("fbar_MW",fbar_MW)
fbar_M31= (D+G)/M31_Totalmass
fbar_M31 = np.round(fbar_M31,3)
print("fbar_M31",fbar_M31)
fbar_M33 = (E)/M33_Totalmass
fbar_M33 = np.round(fbar_M33,3)
print("fbar_M33",fbar_M33)

"""
Calculating mass breakdown of local group by results above
"""
M = C+A+B
print("total halo mass in local group is", np.around(M,3)) # total halo mass in local group
N= F+D+E
print("total disc mass in local group is" , np.around(N,3)) # total disc mass in local group
O = H+G+0
print("total buldge mass in local group is" , np.around(O,3))  # total buldge mass in local group
fbar_Localgroup = np.round((N+O)/M_Localgroup,3)
print("fbar for local group is " , fbar_Localgroup)




    
    

