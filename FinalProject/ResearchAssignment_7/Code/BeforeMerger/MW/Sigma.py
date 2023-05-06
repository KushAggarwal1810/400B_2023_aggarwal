"""

Astr 400 B
Final Project
Kush Aggarwal

"""

"""
The aim of this code to get the output file for each snap id containing the 
velcoity dispersion and black hole mass in x,y,z directions. 

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
from CenterOfMass2 import CenterOfMass          # modified COM class , it now includes volDec.
from buldge_re import eff_buldge
from BM_before import M_sigma



def OrbitCOM(galaxy,start,end,n):
    """function that loops over all the desired snapshots to compute the velcoity dispersion
    and black hole mass in different directions as a function of time.
    
    INPUTS:
        -------
          galaxy: string
          2 word name of the galaxy
          
          start : float
              number of the first snapshot to be read in
          
            end : float
            
            the number of the last snapshot to be read in.
            
            n: float
            
             integer indicating the intervals over which you will return the dispersion 
             velocity and bhm outputs
        
        OUTPUTS: 
        ----------
        saves a text file that has velocity dispersion in different directions
        and black hole mass for each snapshot of time. 
        
    """
    
    # compose the filename for output
    fileout= "%s_" %('sigma')+galaxy+'.txt'  
    
    VolDec =2 
        
    
    # generate the snapshot id sequence 
    
    snap_ids = np.arange(start,end,n)  # giving which snapshots to store
    
    if (snap_ids.size==0):              #checked if the array is empty.
        sys.exit("The array is empty")
        
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    
    Orbit = np.zeros([len(snap_ids),7])
    

    for i, snap_id in enumerate(snap_ids): # loop over 
    
    
        # compose the data filename (be careful about the folder)
        
        ilbl = '000' + str(snap_id)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        filename= "%s_" %(galaxy) + ilbl + '.txt'
        
         
        time1, total1, data1 = Read(filename) 
        t= time1/1000   # converting time grom Myr to Gyr.
        
        e= eff_buldge(filename, VolDec)  # this would give the velocity disperion array
        
        f=  M_sigma(e)                   # this computes the black hole masses in x,y,z
        
        
        Orbit[i]= t.value, *tuple((e)), *tuple((f)) # storing the row at once.
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, Orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'sigx', 'sigy', 'sigz', 'mbx', 'mby', 'mbz'))  

# we identify the merger is happening roughly around 445 snap id.

E= OrbitCOM('MW', 0, 445, 5)