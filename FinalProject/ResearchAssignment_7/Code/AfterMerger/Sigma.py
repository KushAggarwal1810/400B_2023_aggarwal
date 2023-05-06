

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
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterofMass_merger import CenterOfMass          # modified COM class , it now includes volDec.
from buldge_re_m import eff_buldge
from BM_before import M_sigma



def OrbitCOM(galaxy1,galaxy2,start,end,n):
    """function that loops over all the desired snapshots to compute the velcoity dispersion
    and black hole mass in different directions as a function of time for the remnant.
    
    INPUTS:
        -------
          galaxy1: string
          2 word name of the galaxy2
          
          galaxy2: string
          
          2 word name of the galaxy2
          
          
          
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
        and black hole massfor each snapshot of time. 
    """
    
    # compose the filename for output
    fileout= "%s_" %('sigma')+galaxy1+'.txt'  
    
    # fileout = "Orbit_galaxyname.txt"
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    
    VolDec =2 
        
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    
    snap_ids = np.arange(start,end,n)  # giving which snapshots to store
    
    if (snap_ids.size==0):              # BONUS : checked if the array is empty.
        sys.exit("The array is empty")
        
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    
    Orbit = np.zeros([len(snap_ids),7])
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over 
    
    
        # compose the data filename (be careful about the folder)
        
        ilbl = '000' + str(snap_id)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        filename1= "%s_" %(galaxy1) + ilbl + '.txt'
        filename2= "%s_" %(galaxy2) + ilbl + '.txt'
        
        e= eff_buldge(filename1,filename2,VolDec)
        
        time1, total1, data1 = Read(filename1) 
        t= time1/1000   # converting time grom Myr to Gyr.
        
        e= eff_buldge(filename1,filename2,VolDec)
        f=  M_sigma(e)
        
        
        Orbit[i]= t.value, *tuple((e)), *tuple((f)) # storing the row at once.
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, Orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'sigx', 'sigy', 'sigz', 'mbx', 'mby', 'mbz'))  



E= OrbitCOM('MW','M31', 445, 801, 5)