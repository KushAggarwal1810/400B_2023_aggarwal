# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:36:06 2023

@author: aggar
"""

import numpy as np
import astropy.units as u

# This function opens and reads an input file
# inputs : filename ( preferably in .txt format)
# returns : the time , total number of particles , particle type, mass, x,y,z, vx,vy,vz columns as a data array

def Read(filename):
    file = open(filename,'r' )      # opens the input text file
    line1 = file.readline()         # reads the first line of the text file
    label, value = line1.split()    # Stores the label and value as different things
    time = float(value)*u.Myr       # converts the 'value' into float and assigns Myr units to it.
    line2 = file.readline()         # reads the second line of the text file
    label, value= line2.split()     # Stores the label and value as different things
    total_particles = float(value)  # converts the 'value' into float  
    file.close()                    # closes the file
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3) 
    
    #this reads remainder of the file and stores each value under the repsective
    #label in the data array
    
    
    # print(data['type'][1])  # this is the example used to test the function
    
    return time,total_particles,data  # return statement that returns time, total particles and rest of the parameters as a data array








