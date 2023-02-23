

# Homework 6 Template
# G. Besla & R. Li




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
from CenterOfMass2 import CenterOfMass          # modified COM class , it now includes volDec.




def OrbitCOM(galaxy,start,end,n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
          galaxy: string
          2 word name of the galaxy
          
          start : float
              number of the first snapshot to be read in
          
            end : float
            
            the number of the last snapshot to be read in.
            
            n: float
            
             integer indicating the intervals over which you will return the COM
    outputs: 
        
        saves a text file that has com postions and velocities
        for each snapshot of time. 
    """
    
    # compose the filename for output
    fileout= "%s_" %('Orbit')+galaxy+'.txt'  
    
    # fileout = "Orbit_galaxyname.txt"
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    if ( galaxy== 'M33'):    # here I am distinguishing M33 from other galaxies
        delta = 0.1
        VolDec =4
    
    else :
        delta = 0.1
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
        filename= "%s_" %(galaxy) + ilbl + '.txt'
        
        # Initialize an instance of CenterOfMass class, using disk particles ie 2
        
        COM = CenterOfMass(filename, 2)

        # Store the COM pos and vel. Remember that now COM_P required VolDec
         
        COM_p = COM.COM_P(0.1,VolDec)   # using our COM_P function from COM class to compute COM position coord.
        COM_v = COM.COM_V(COM_p[0], COM_p[1], COM_p[2]) # calculating COM_V velcoties coordinates.
        
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        
        t= COM.time/1000   # converting time grom Myr to Gyr.
        
        Orbit[i]= t.value, *tuple((COM_p).value), *tuple((COM_v).value) # storing the row at once.
        
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

# E= OrbitCOM('MW', 0, 801, 5)


# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt

"""
Now, we have run the above cell, and have got the orbit files for each galaxy, so we will 
read the files for data now...
"""
file = open('Orbit_MW.txt','r' )                
data = np.genfromtxt('Orbit_MW.txt',dtype=None,names=True)  # storing all labels to data

t_MW = data['t']
x_MW = data['x']
y_MW= data['y']
z_MW= data['z']
vx_MW= data['vx']
vy_MW= data['vy']
vz_MW= data['vz']

file = open('Orbit_M31.txt','r' )
data = np.genfromtxt('Orbit_M31.txt',dtype=None,names=True) 

t_M31 = data['t']
x_M31 = data['x']
y_M31= data['y']
z_M31= data['z']
vx_M31= data['vx']
vy_M31= data['vy']
vz_M31= data['vz']

file = open('Orbit_M33.txt','r' )
data = np.genfromtxt('Orbit_M33.txt',dtype=None,names=True)
t_M33 = data['t']
x_M33 = data['x']
y_M33= data['y']
z_M33= data['z']
vx_M33= data['vx']
vy_M33= data['vy']
vz_M33= data['vz'] 

# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  


def mag(a,b):
    diff = np.array([])
    for i in range(len(a)):
        t= abs(a[i]-b[i])
        diff = np.append(diff,t)
    v=LA.norm(diff)
    return v   
        



# Determine the magnitude of the relative position and velocities 

# of MW and M31

rel_sep1= np.array([])
for i in range(len(t_MW)):
    a = np.array([x_MW[i],y_MW[i],z_MW[i]])
    b= np.array([x_M31[i],y_M31[i],z_M31[i]])
    c= mag(a,b)
    rel_sep1= np.append(rel_sep1,c)
    

# P= min(rel_sep1)   # finding where MW M31 merge
# index = np.where(rel_sep1==P)

# print(t_MW[90])
         
rel_sep2= np.array([])
for i in range(len(t_MW)):
    a = np.array([vx_MW[i],vy_MW[i],vz_MW[i]])
    b= np.array([vx_M31[i],vy_M31[i],vz_M31[i]])
    c= mag(a,b)
    rel_sep2= np.append(rel_sep2,c)   

# of M33 and M31

rel_sep3= np.array([])
for i in range(len(t_MW)):
    a = np.array([x_M33[i],y_M33[i],z_M33[i]])
    b= np.array([x_M31[i],y_M31[i],z_M31[i]])
    c= mag(a,b)
    rel_sep3= np.append(rel_sep3,c)

rel_sep4= np.array([])
for i in range(len(t_MW)):
    a = np.array([vx_M33[i],vy_M33[i],vz_M33[i]])
    b= np.array([vx_M31[i],vy_M31[i],vz_M31[i]])
    c= mag(a,b)
    rel_sep4= np.append(rel_sep4,c)   


# Plot the Orbit of the galaxies 
#################################

"""
Plotting Relative Separation between MW and M31
"""
# plt.semilogy(t_MW, rel_sep1,color='orange',label=r'MW-M33')
plt.plot(t_MW, rel_sep1,color='orange',label=r'MW-M33')
plt.xlabel('Time(Gyr)')
plt.ylabel('Separation(kpc)')
plt.title('Figure for Relative Separation between MW and M31')
plt.show()

"""
Plotting Relative Separation between M33 and M31
"""

plt.plot(t_MW, rel_sep3,color='red',label=r'M33-M31')
plt.xlabel('Time(Gyr)')
plt.ylabel('Separation(kpc)')
plt.title('Figure for Relative Separation between M33 and M31')
plt.show()




# Plot the orbital velocities of the galaxies 
#################################

"""
Plotting Relative velocity Between between MW and M31
"""

plt.plot(t_MW, rel_sep2,color='green',label=r'MW-M31')
plt.xlabel('Time(Gyr)')
plt.ylabel('Separation(kpc)')
plt.title('Figure for Relative velocity between MW and M31')
plt.show()

"""
Plotting Relative velocity Between between M31 and M33
"""

plt.plot(t_MW, rel_sep4,color='blue',label=r'M31-M33')
plt.xlabel('Time(Gyr)')
plt.ylabel('Separation(kpc)')
plt.title('Figure for Relative velocity between M31 and M33')
plt.show()
