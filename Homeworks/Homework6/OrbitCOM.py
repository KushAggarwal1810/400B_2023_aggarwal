

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
import sys

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass




def OrbitCOM(galaxy,start,end,n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
          
    outputs: 
    """
    
    # compose the filename for output
    fileout= "%s_" %('Orbit')+galaxy+'.txt'
    
    # fileout = "Orbit_galaxyname.txt"
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    if ( galaxy== 'M33'):
        delta = 0.1
        VolDec =4
    
    else :
        delta = 0.1
        VolDec =2 
        
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    
    snap_ids = np.arange(start,end,n)
    
    if (snap_ids.size==0):
        sys.exit("The array is empty")
        
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    
    Orbit = np.zeros([len(snap_ids),7])
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over 
    
    
        # compose the data filename (be careful about the folder)
        
        ilbl = '000' + str(snap_id)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        filename= "%s_" %(galaxy) + ilbl + '.txt'
        
        # Initialize an instance of CenterOfMass class, using disk particles
        
        COM = CenterOfMass(filename, 2)

        # Store the COM pos and vel. Remember that now COM_P required VolDec
        
        COM_p = COM.COM_P(0.1,VolDec)
        COM_v = COM.COM_V(COM_p[0], COM_p[1], COM_p[2])
        
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        
        t= COM.time/1000
        
        Orbit[i]= t, *tuple(COM_p.value), *tuple(COM_v.value)
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, Orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))




# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 




# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt




# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  




# Determine the magnitude of the relative position and velocities 

# of MW and M31

# of M33 and M31




# Plot the Orbit of the galaxies 
#################################




# Plot the orbital velocities of the galaxies 
#################################

