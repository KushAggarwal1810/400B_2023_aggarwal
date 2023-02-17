# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:48:23 2023

@author: aggar
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 13:25:37 2023

@author: aggar
"""

import numpy as np
import astropy.units as u
import astropy.table as tbl
from numpy import linalg as LA
import matplotlib.pyplot as plt
from scipy.integrate import quad # For integration
from ReadFile import Read
from CenterOfMass import  CenterOfMass

class MassProfile:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot

    def __init__(self, galaxy, Snap):
        
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            galaxy: A String with Galaxy Name, e.g. “MW”, “M31” or “M33”
            
            snap: Snapshot number, e.g. 0, 1, etc
        '''
        
        # add a string of the filenumber to the value “000”
        ilbl = '000' + str(Snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename= "%s_" %(galaxy) + ilbl + '.txt'
        
        # read data in the given file using Read
        self.time, self.total, self.data = Read(self.filename) 
  
       
        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m']
        # write your own code to complete this for positions and velocities
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        self.vx = self.data['vx']
        self.vy = self.data['vy']
        self.vz = self.data['vz']
        # print(self.x[1])
        
        self.gname= galaxy
        
        # print(self.gname)

    def MassEnclosed(self,ptype,r):
        """
        
        Method to compute the mass enclosed in a specific component of a galaxy
        
        inputs:
           ptype : `int; 1, 2, or 3`
               particle type to use for COM calculations
        
           r : array of radii
           
         Ouput :
             
             m_enc : astropy qty
                mass enclosed of a particular component at diff. r values
        """
        MW_COM = CenterOfMass(self.filename, ptype)
        MW_COM_p_MW = MW_COM.COM_P(0.1)
        
        # print(MW_COM_p_MW)
        
        #create an array to store indexes of particles of desired Ptype                                
        index = np.where(self.data['type'] == ptype)
        
        m1 = self.m[index]
        x1= self.x[index]
        y1 = self.y[index]
        z1= self.z[index]
        
        # print(z1)
        
        xM =  abs(MW_COM_p_MW[0]-x1)  # converting coordinates to COM frame
        yM =  abs(MW_COM_p_MW[1]-y1)
        zM = abs(MW_COM_p_MW[2]-z1)
        
        rM = np.sqrt(xM**2  + yM**2 + zM**2) # this is the radial vector of ptype particles wrt COM ie the origin
        
        m_enc = np.array([])
        
        
        for i in range(len(r)):
            
            index2 = np.where(rM.value <= r[i].value)
            t= m1[index2]
            t_sum = np.sum(t)*1e10               # converting to 1M_sun 
            
            m_enc = np.append(m_enc, t_sum)
            
        m_enc = m_enc * u.M_sun
            
        
        return m_enc
            
    def MassEnclosedTotal(self,r):
        
        
        """
        
        Method to compute total enclosed mass in a galaxy
        
        input parameters :
            
            r = array of radii
            
        Ouput 
        
        m_total : astropy qty
            total mass enlcosed in the galaxy ( halo + disk + buldge) at diff r values
        """
        
        if (self.gname=='M33'):    # using if / else, to incorporate the caveat that M33 doesn't have a buldge.
            
            # buldge_mass = 0
            
            disk_mass = self.MassEnclosed(2,r)
            
            halo_mass = self.MassEnclosed(1,r)
            
            m_total= np.add(halo_mass,disk_mass)
        
            
        else :
            
            halo_mass = self.MassEnclosed(1,r)
            
            disk_mass = self.MassEnclosed(2,r)
            
            buldge_mass = self.MassEnclosed(3,r)
            
       
            
            m_total = np.add(halo_mass,disk_mass,buldge_mass)
            
            m_total = m_total 
            
        
        return m_total
         
    def HernquistMass(self,r,a,Mhalo):
        
        """
        
        Method to fit the hernquist mass profile 
        
        Inputs :
            
            r : array of radius in kpc
            a : scale length 
            
            Mhalo : total halo mass of the gaalxy in solar mass units
            
        Ouputs : astropy qty
        
            Total halo mass enclosed in the gaalxy as per hernquist profile at diff r values in solar mass units
        
        """
        
          
        m_enc_h = np.array([])
       
        
               
        for i in range(len(r)):
            
            t =  (Mhalo.value*r[i].value**2)/(a.value+r[i].value)**2   # using the hernquit mass profile formula
            
           
            m_enc_h= np.append(m_enc_h,t)
            
        m_enc_h=  m_enc_h*u.M_sun
        
        # print(m_enc_h)
            
                
        return m_enc_h
    
    
            
    def CircularVelocity(self,ptype,r):
        
        """
        
        Method to compute the Circular Velocity of the particular galaxy component
        
        Inputs :
            
            r : array of radii in kpc
            ptype : particle type : 1= halo , 2 = disk , 3= buldge
        
        Outputs :
            
            v_cir : astropy qty
            
            circular velocity of the inout ptype at diferent radii values in Km/s
        """
        
        from astropy.constants import G
        
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        m =  self. MassEnclosed(ptype,r)
        v_cir = np.sqrt(G*m/r)
        
        return v_cir
        
    
    def TotalCircularVelocity(self,r):
        
        """
        
        Method to compute the total Circular Velocity at each radii due to all galaxy components.
        
        Inputs :
            
            r : array of radii in kpc
          
        
        Outputs :
            
            v_cir_total : astropy qty
            
            total circular velocity of  at diferent radii values in km/s
        """
        
        from astropy.constants import G
        
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
    
        
        m_total =  self.MassEnclosedTotal(r)  # we cannot just do sum of imdvidual vecloities as their is a sqrt
        
        v_cir_total = np.sqrt(G*m_total/r)
        
        # print(v_cir_total)
        
        return v_cir_total
    
    def HernquistVCirc(self,r,a,Mhalo):
        
        """
        Method to compute the Hernquist circular velocity , which is the circular speed computed 
        using hernquist mass profile
        
        Inputs :
            
            r = radii in kpc
            
            a=  scale radius
            
            Mhalo : total halo mass of thr galaxy in solar mass units
            
        Outputs :
            
            v_cir_h : astropy qty
            
            Hernquist circular speed in km/s computed using hernquist mass profile at diff radii
        
        
        
        """
        from astropy.constants import G
        
        r = r
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        m =  self.HernquistMass(r,a, Mhalo)
        v_cir_h = np.sqrt(G.value*m/r)
        
        return v_cir_h

        
    def total_halomass(self):    
        
        """
        Method for computing total halo mass enclosed in the galaxy 
        
        inputs : none, just the initialized qty in class
        
        outputs :
            
            mt = astropy qty 
            
            total halo mass of the galaxy
        """
    
     
        index = np.where(self.data['type']==1)
        mass = self.m[index]
        
        mt = np.sum(mass)*1e10
        
        mt = mt * u.M_sun
        
        return mt
        
        
        


if __name__ == '__main__' : 

    """ Mass profiles for MW"""

    MW = MassProfile("MW", 0)  
    
    r = np.arange(0.1, 30, 1) 

    r = r*u.kpc    
    
    Mass_Encl_halo = MW.MassEnclosed(1, r)
    
    
    Mass_Encl_disk = MW.MassEnclosed(2, r)
     
    Mass_Encl_buldge=  MW.MassEnclosed(3, r)
    
  
    
    Mass_Enc_Total= MW.MassEnclosedTotal(r) 
    
    M_halo = MW.total_halomass()
    
    a= 62.5*u.kpc
    
    Mass_Her = MW.HernquistMass(r, a,M_halo )
    
  
    
    plt.plot(r,Mass_Encl_buldge,color='orange', linestyle='-', linewidth=3, 
                label=r'Buldge Mass')
  
    
    plt.plot(r,Mass_Encl_halo,color='red', linestyle=':', linewidth=3, 
                label=r'halo mass')
    
    plt.plot(r,Mass_Encl_disk,color='purple', linestyle='solid', linewidth=3, 
                label=r'disk mass')
    
    
    plt.plot(r,Mass_Enc_Total,color='green', linestyle='-.', linewidth=3, 
                label=r'total enclosed mass')
    
    plt.plot(r,Mass_Her,color='blue', linestyle='--', linewidth=3, 
                label=r'Herniquist mass profile ; scale length = $62.5$')
    
    plt.title(" Mass profile for MW")
    plt.xlabel("radius in kpc")
    plt.ylabel (r'Mass profiles in solar mass units')
    
    plt.yscale('log')
    plt.legend()
    plt.show()
    
    
    """ Mass profiles for M31"""
    
    M31 = MassProfile("M31", 0)  
    
    r = np.arange(0.1, 30, 1) 

    r = r* u.kpc    
    
    Mass_Encl_halo = M31.MassEnclosed(1, r)
    
    
    Mass_Encl_disk = M31.MassEnclosed(2, r)
     
    Mass_Encl_buldge=  M31.MassEnclosed(3, r)
    
    
    
    Mass_Enc_Total= M31.MassEnclosedTotal(r) 
    
    M_halo = M31.total_halomass()
    
    a= 61.5* u.kpc
    
    Mass_Her = M31.HernquistMass(r, a,M_halo )
    
    
    
    plt.plot(r,Mass_Encl_buldge,color='orange', linestyle='-', linewidth=3, 
                label=r'Buldge Mass')
    
    
    plt.plot(r,Mass_Encl_halo,color='red', linestyle=':', linewidth=3, 
                label=r'halo mass')
    
    plt.plot(r,Mass_Encl_disk,color='purple', linestyle='solid', linewidth=3, 
                label=r'disk mass')
    
    
    plt.plot(r,Mass_Enc_Total,color='green', linestyle='-.', linewidth=3, 
                label=r'total enclosed mass')
    
    plt.plot(r,Mass_Her,color='blue', linestyle='--', linewidth=3, 
                label=r'Herniquist mass profile ; scale length = $61.5$')
    
     
    plt.title(" Mass profile for M31")
    plt.xlabel("radius in kpc")
    plt.ylabel (r'Mass profiles in solar mass units')
    
    plt.yscale('log')
    plt.legend()
    plt.show()
    
 
     
    

    """ Mass profiles for M33"""
    
    M33 = MassProfile("M33", 0)  
    
    r = np.arange(0.1, 30, 1) 

    r= r*u.kpc    
    
    Mass_Encl_halo = M33.MassEnclosed(1, r)
    
    
    Mass_Encl_disk = M33.MassEnclosed(2, r)
     
    

    Mass_Enc_Total= M33.MassEnclosedTotal(r) 
    
    M_halo = M33.total_halomass()
    
    a= 26* u.kpc
    
    Mass_Her = M33.HernquistMass(r, a,M_halo )
    
    
    
   
    
    
    plt.plot(r,Mass_Encl_halo,color='red', linestyle=':', linewidth=3, 
                label=r'halo mass')
    
    plt.plot(r,Mass_Encl_disk,color='purple', linestyle='solid', linewidth=3, 
                label=r'disk mass')
    
    
    plt.plot(r,Mass_Enc_Total,color='green', linestyle='-.', linewidth=3, 
                label=r'total enclosed mass')
    
    plt.plot(r,Mass_Her,color='blue', linestyle='--', linewidth=3, 
                label=r'Herniquist mass profile ; scale length = $26$')
    
     
    plt.title(" Mass profile for M33")
    plt.xlabel("radius in kpc")
    plt.ylabel (r'Mass profiles in solar mass units')
    
    plt.yscale('log')
    plt.legend()
    plt.show()
    
       
    """ Circular velocity profiles for MW"""
    
    MW = MassProfile("MW", 0)  
    
    r = np.arange(0.1, 30, 1) 

    r= r*u.kpc    
    
    Circular_halo = MW.CircularVelocity(1, r)
    
    
    Circular_disk = MW.CircularVelocity(2, r)
     
    Circular_buldge = MW.CircularVelocity(3, r)
    
    
    
    Circular_Total= MW.TotalCircularVelocity(r)
    
    M_halo = MW.total_halomass()
    
    a= 62*u.kpc
    
    Cir_Her = MW.HernquistVCirc(r, a,M_halo )
    
    
    
    plt.plot(r,Circular_buldge,color='orange', linestyle='-', linewidth=3, 
                label=r'Buldge Mass')
    
    
    plt.plot(r, Circular_halo,color='red', linestyle=':', linewidth=3, 
                label=r'halo mass')
    
    plt.plot(r, Circular_disk,color='purple', linestyle='solid', linewidth=3, 
                label=r'disk mass')
    
    
    plt.plot(r, Circular_Total,color='green', linestyle='-.', linewidth=3, 
                label=r'total enclosed mass')
    
    plt.plot(r, Cir_Her,color='blue', linestyle='--', linewidth=3, 
                label=r'Herniquist mass profile ; scale radius = $62$')
    
    plt.title(" Circular Velocity profile for MW")
    plt.xlabel("radius in kpc")
    plt.ylabel (r'Circular velocity profiles in $\frac{Km}{s}$')
    
    
    plt.yscale('log')
    plt.legend(loc = 'lower center', fontsize ='9')
    plt.show()
    
       
    """ Circular velocity profiles for M31"""
    
    M31 = MassProfile("M31", 0)  
    
    r = np.arange(0.1, 30, 1)     
    
    r=  r*u.kpc
    
    Circular_halo = M31.CircularVelocity(1, r)
    
    
    Circular_disk = M31.CircularVelocity(2, r)
     
    Circular_buldge = M31.CircularVelocity(3, r)
    
    
    
    Circular_Total= M31.TotalCircularVelocity(r)
    
    M_halo = M31.total_halomass()
    
    a= 62 * u.kpc
    
    Cir_Her = M31.HernquistVCirc(r, a,M_halo )
    
    
    
    plt.plot(r,Circular_buldge,color='orange', linestyle='-', linewidth=3, 
                label=r'Buldge Mass')
    
    
    plt.plot(r, Circular_halo,color='red', linestyle=':', linewidth=3, 
                label=r'halo mass')
    
    plt.plot(r, Circular_disk,color='purple', linestyle='solid', linewidth=3, 
                label=r'disk mass')
    
    
    plt.plot(r, Circular_Total,color='green', linestyle='-.', linewidth=3, 
                label=r'total enclosed mass')
    
    plt.plot(r, Cir_Her,color='blue', linestyle='--', linewidth=3, 
                label=r'Herniquist mass profile ; scale radius = $62$')
    
    plt.title(" Circular Velocity profile for M31")
    plt.xlabel("radius in kpc")
    plt.ylabel (r'Circular velocity profiles in $\frac{Km}{s}$')
    
    
    plt.yscale('log')
    plt.legend()
    plt.show()
    
    
     
    
    
    """ Circular velocity profiles for M33"""
    
    M33 = MassProfile("M33", 0)  
    
    r = np.arange(0.1, 30, 1)    
    
    r=  r*u.kpc
    
    Circular_halo = M33.CircularVelocity(1, r)
    
    
    Circular_disk = M33.CircularVelocity(2, r)
     

    
    
    
    Circular_Total= M33.TotalCircularVelocity(r)
    
    M_halo = M33.total_halomass()
    
    a= 26.5* u.kpc
    
    Cir_Her = M33.HernquistVCirc(r, a,M_halo )
    
    
    plt.plot(r, Circular_halo,color='red', linestyle=':', linewidth=3, 
                label=r'halo mass')
    
    plt.plot(r, Circular_disk,color='purple', linestyle='solid', linewidth=3, 
                label=r'disk mass')
    
    
    plt.plot(r, Circular_Total,color='green', linestyle='-.', linewidth=3, 
                label=r'total enclosed mass')
    
    plt.plot(r, Cir_Her,color='blue', linestyle='--', linewidth=3, 
                label=r'Herniquist mass profile ; scale radius = $ 26.5$')
    
     
    plt.title(" Circular Velocity profile for M33")
    plt.xlabel("radius in kpc")
    plt.ylabel (r'Circular velocity profiles in $\frac{Km}{s}$')
    
    
    plt.yscale('log')
    plt.legend()
    plt.show()
    
    
    
    
    
     
    
    
    
    
    
     