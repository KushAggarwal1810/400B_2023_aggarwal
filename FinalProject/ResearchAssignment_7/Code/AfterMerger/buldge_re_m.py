
"""

Astr 400 B
Final Project
Kush Aggarwal

"""

"""
The aim of this code to find the effective buldge radii of the input galaxy file
and map the correpsonding buldge particles, then it computes the velcoity dispersions
in x,y and z directions.  This code is different from before merger ans afte rmerger we need 
to consider the remanant so we need to redefine our effective buldge radii to
include the buldge particles of both MW and M31


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
# my modules
from ReadFile import Read 
from CenterofMass_merger import CenterOfMass 


def eff_buldge(filename1,filename2,VolDec):
    
        """
        This function calculates the effective buldge radii for the remnanat and returns
        the velocity dispersion in x,y,z directions
        
        INPUTS:
        -------
        
        Filename 1: filename in .txt of 1st galaxy
        
        Filename 2: filename in .txt of 2nd galaxy
        
        Voldec : The voldec of the two galaxies ( MW, M31) is same here, so we define
        a common voldec
        
        OUTPUTS:
            
            sig : array of velocity disperisons in x,y,z directions.
            
        """
        
        time1, total1, data1 = Read(filename1) 
        
        """
        Storing buldge particles data from first file
        """
        index1 = np.where(data1['type'] == 3)
       
      
       
        m1 = data1['m'][index1]
        x1 = data1['x'][index1]*u.kpc
        y1 = data1['y'][index1]*u.kpc
        z1 = data1['z'][index1]*u.kpc
        vx1 = data1['vx'][index1]
        vy1 = data1['vy'][index1]
        vz1 = data1['vz'][index1]
        
        time2, total2, data2 = Read(filename2) 
        
        """
        Storing buldge particles data from second file
        """
        
        index2 = np.where(data2['type'] == 3)
        
        m2 = data2['m'][index2]
        x2=  data2['x'][index2]*u.kpc
        y2 = data2['y'][index2]*u.kpc
        z2=  data2['z'][index2]*u.kpc
        vx2 = data2['vx'][index2]
        vy2 = data2['vy'][index2]
        vz2 = data2['vz'][index2]
        
        """
        We use np.concatenate to get the total buldge particles of MW, M31
        """
        
        m= np.concatenate((m1,m2))
        x= np.concatenate((x1,x2))
        y= np.concatenate((y1,y2))
        z= np.concatenate((z1,z2))
        vx= np.concatenate((vx1,vx2))
        vy= np.concatenate((vy1,vy2))
        vz= np.concatenate((vz1,vz2))
        
        COM = CenterOfMass(filename1,filename2, 2)     # finding center of mass of remnant wrt disk
        COM_p = COM.COM_P(0.1,VolDec)
        
        xM =  abs(COM_p[0]-x)  # converting coordinates to COM frame
        yM =  abs(COM_p[1]-y)
        zM = abs(COM_p[2]-z)
        
        r = np.arange(0.1, 30, 1)    # radii array
        
        rM = np.sqrt(xM**2  + yM**2 + zM**2)
        
        m_enc = np.array([])
        
        """
        
        We need to find the total buldge mass enclosed
        """
        
        for i in range(len(r)):
            
            index4 = np.where(rM.value <= r[i])  # rM is an array
            t= m[index4]
            t_sum = np.sum(t)*1e10               # converting to 1M_sun 
            
            m_enc = np.append(m_enc, t_sum)
            
        m_enc = m_enc * u.M_sun
        
        t= max(m_enc)   # this is total buldge mass enclosed
        
        y5=[]
        
        
        """
        Now, we will find the index that contains half the buldge particle mass
        
        """
        for i in range(len(r)):
            
            if (m_enc[i]>= t/2):
                y5.append(i)

        g= y5[0]     # required index
    
        r_d= r[g]     # eff radii : radii that contains half the buldge mass
            
        
        """
        
        Storing the buldge particles in eff radii
        """
        index3 = np.where(rM.value <= r_d)
        
        # print(index3)
        
        m3 = m[index3] 
        # print(m3) # these are the parameters of buldge particles inside the r_eff
        x3=  x[index3]
        y3 = y[index3]
        z3=  z[index3]
        vx3 = vx[index3]
        vy3 = vy[index3]
        vz3 = vz[index3]
        
        # print(len(vx1))
        # print(len(vx2))
        """
        We know compute the velocity dispersion using inbuilt np.std function
        """
        sigma_x= np.std(vx3)
        sigma_y= np.std(vy3)
        sigma_z= np.std(vz3)
        
        sig= np.array([sigma_x,sigma_y,sigma_z])
        
        return sig
        
            
            
     


