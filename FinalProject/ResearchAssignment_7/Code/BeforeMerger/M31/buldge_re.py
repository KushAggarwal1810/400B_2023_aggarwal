"""

Astr 400 B
Final Project
Kush Aggarwal

"""

"""
The aim of this code to find the effective buldge radii of the input galaxy file
and map the correpsonding buldge particles, then it computes the velcoity dispersions
in x,y and z directions. 

"""
# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
import sys
from numpy import linalg as LA
# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
from ReadFile import Read
from CenterOfMass2 import CenterOfMass   


def eff_buldge(filename,VolDec):
    
        """
        Inputs:
        -------
        
        filename : snap id filename in .txt format
        Voldec   : Voldec of the Galaxy under consideration
        
        OUTPUTS:
        -------
        
        sig : array with x dispersion vel, y disp vel and z disp vel.
        
        
        """
        
        time, total, data = Read(filename) 
        
        """
        We read in each of the columns
        """
        m = data['m']
        x = data['x']*u.kpc
        y = data['y']*u.kpc
        z = data['z']*u.kpc
        vx = data['vx']
        vy = data['vy']
        vz = data['vz']
     
        COM = CenterOfMass(filename, 2)        # Computes Center of Mass with respect to disk
        COM_p = COM.COM_P(0.1,VolDec)
        index = np.where(data['type'] == 3)    # Index Buldge Particles
        
        """
        Now we store only the buldge particles
        """
        
        m1 = m[index]
        x1=  x[index]
        y1 = y[index]
        z1=  z[index]
        vx1 = vx[index]
        vy1 = vy[index]
        vz1 = vz[index]
        
        # print(z1)
        
        xM =  abs(COM_p[0]-x1)  # converting coordinates to COM frame
        yM =  abs(COM_p[1]-y1)
        zM = abs(COM_p[2]-z1)
        
        r = np.arange(0.1, 30, 1)    # radii array
        
        rM = np.sqrt(xM**2  + yM**2 + zM**2)    # distance from center of mass
        
        m_enc = np.array([])
        
        """
        We know define the enclosed mass in order to get the total mass of buldge
        """
        for i in range(len(r)):
            
            index2 = np.where(rM.value <= r[i])  # rM is an array
            t= m1[index2]
            t_sum = np.sum(t)*1e10               # converting to 1M_sun 
            
            m_enc = np.append(m_enc, t_sum)
            
        m_enc = m_enc * u.M_sun
        
        t= max(m_enc)                           # This is the total mass                
        
        y5=[]
        
        for i in range(len(r)):
            
            if (m_enc[i]>= t/2):               # Finding the index where we have half the enclosed mass
                y5.append(i)

        g= y5[0]                               # This is the required index
    
        r_d= r[g]                              #  This is the corresponding radii in the radius array, uptil which we have have the enclosed mass
            
        index2 = np.where(rM.value <= r_d)
        
        m2 = m1[index2]      # these are the parameters of buldge particles inside the r_eff
        x2=  x1[index2]
        y2 = y1[index2]
        z2=  z1[index2]
        vx2 = vx1[index2]
        vy2 = vy1[index2]
        vz2 = vz1[index2]
        
        # print(len(vx1))
        # print(len(vx2))
        
        """
        We know compute the velocity dispersion using inbuilt np.std function
        """
        sigma_x= np.std(vx2)
        sigma_y= np.std(vy2)
        sigma_z= np.std(vz2)
        
        sig= np.array([sigma_x,sigma_y,sigma_z])
        
        return sig
        
            
            
       




