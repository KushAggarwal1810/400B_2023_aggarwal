"""

Astr 400 B
Final Project
Kush Aggarwal

"""

"""
The aim of this code to modify the center of mass code made in class to be able
to calculate the combined center of mass of two galaxies.

"""

# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl
from numpy import linalg as LA
from ReadFile import Read




class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot

    def __init__(self, filename1,filename2, ptype):
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            filename1 : `str`
                snapshot file of galaxy 1
            
            filename1 : `str`
                snapshot file of galaxy 2
            ptype : `int; 1, 2, or 3`
                particle type to use for COM calculations
        '''
     
        """
        We read in data from 1st file
        """
        # read data in the given file using Read
        self.time1, self.total1, self.data1 = Read(filename1)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index1 = np.where(self.data1['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m1 = self.data1['m'][self.index1]
        # write your own code to complete this for positions and velocities
        self.x1 = self.data1['x'][self.index1]
        self.y1 = self.data1['y'][self.index1]
        self.z1 = self.data1['z'][self.index1]
        self.vx1 = self.data1['vx'][self.index1]
        self.vy1 = self.data1['vy'][self.index1]
        self.vz1 = self.data1['vz'][self.index1]
        
        
        """
        We read in data from 2nd file
        """
        # read data in the given file using Read
        self.time2, self.total2, self.data2 = Read(filename2)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index2 = np.where(self.data2['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m2 = self.data2['m'][self.index2]
        # write your own code to complete this for positions and velocities
        self.x2 = self.data2['x'][self.index2]
        self.y2 = self.data2['y'][self.index2]
        self.z2 = self.data2['z'][self.index2]
        self.vx2 = self.data2['vx'][self.index2]
        self.vy2 = self.data2['vy'][self.index2]
        self.vz2 = self.data2['vz'][self.index2]


        """
        We combine the two data
        """
        self.m= np.concatenate((self.m1,self.m2))
        self.x= np.concatenate((self.x1,self.x2))
        self.y= np.concatenate((self.y1,self.y2))
        self.z= np.concatenate((self.z1,self.z2))
        self.vx= np.concatenate((self.vx1,self.vx2))
        self.vy= np.concatenate((self.vy1,self.vy2))
        self.vz= np.concatenate((self.vz1,self.vz2))
        
        
        """
        Rest of the code is same as class
        """


    def COMdefine(self,a,b,c,m):
        ''' Method to compute the COM of a generic vector quantity by direct weighted averaging.
        
        PARAMETERS
        ----------
        a : `float or np.ndarray of floats`
            first vector component
        b : `float or np.ndarray of floats`
            second vector component
        c : `float or np.ndarray of floats`
            third vector component
        m : `float or np.ndarray of floats`
            particle masses
        
        RETURNS
        -------
        a_com : `float`
            first component on the COM vector
        b_com : `float`
            second component on the COM vector
        c_com : `float`
            third component on the COM vector
        '''
        # write your own code to compute the generic COM 
        #using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        a_com = np.dot(m,a)/np.sum(m)
        # ycomponent Center of mass
        b_com = np.dot(m,b)/np.sum(m)
        # zcomponent Center of mass
        c_com = np.dot(m,c)/np.sum(m)
        
        
        
        # return the 3 components separately
        return a_com,b_com,c_com
    
    
    def COM_P(self, delta, volDec):
        '''Method to compute the position of the center of mass of the galaxy 
        using the shrinking-sphere method.
        PARAMETERS
        ----------
        delta : `float, optional`
            error tolerance in kpc. Default is 0.1 kpc
        
        RETURNS
        ----------
        p_COM : `np.ndarray of astropy.Quantity'
            3-D position of the center of mass in kpc
        '''                                                                     

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        # write your own code below
        r_COM = np.sqrt(x_COM**2 + y_COM**2 + z_COM**2)


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        x_new = self.x -x_COM
        y_new = self.y -y_COM
        z_new = self.z -z_COM
        r_new = np.sqrt(x_new**2 + y_new**2 +z_new**2)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        r_max = max(r_new)/volDec
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        change = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (change > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(r_new <= r_max)
            x3 =  self.x[index2]
            y3 =  self.y[index2]
            z3 =  self.z[index2]
            m3 =  self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x3, y3,z3,m3)
            # compute the new 3D COM position
            # write your own code below
            r_COM2 = np.sqrt(x_COM2**2 + y_COM2**2 + z_COM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            change = np.abs(r_COM - r_COM2)
            # uncomment the following line if you want to check this                                                                                               
            # print ("CHANGE = ", change)                                                                                     

            # Before loop continues, reset : r_max, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            r_max /= volDec
            # check this.                                                                                              
            # print ("maxR", r_max)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            x_new = self.x -x_COM2
            y_new = self.y -y_COM2
            z_new = self.z -z_COM2
            r_new = np.sqrt(x_new**2 + y_new**2 +z_new**2) 

            # set the center of mass positions to the refined values                                                   
            x_COM = x_COM2
            y_COM = y_COM2
            z_COM = z_COM2
            r_COM = r_COM2

            # create an array (np.array) to store the COM position                                                                                                                                                       
        p_COM = np.array([x_COM, y_COM, z_COM])
            
            
                
        # set the correct units using astropy and round all values
        # and then return the COM positon vector
        # write your own code below
        p_COM = np.round(p_COM,2)
        p_COM = p_COM*u.kpc              # converting to kpc
             
        return p_COM
        
        
    def COM_V(self, x_COM, y_COM, z_COM):
        ''' Method to compute the center of mass velocity based on the center of mass
        position.
        PARAMETERS
        ----------
        x_COM : 'astropy quantity'
            The x component of the center of mass in kpc
        y_COM : 'astropy quantity'
            The y component of the center of mass in kpc
        z_COM : 'astropy quantity'
            The z component of the center of mass in kpc
            
        RETURNS
        -------
        v_COM : `np.ndarray of astropy.Quantity'
            3-D velocity of the center of mass in km/s
        '''
        
        # the max distance from the center that we will use to determine 
        #the center of mass velocity                   
        rv_max = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position (x_COM, y_COM, z_COM)
        # write your own code below
        xV = self.x -x_COM.value    # as self.x is unitless , so only need value of xcom, ie make it unitless.
        yV = self.y -y_COM.value
        zV = self.z -z_COM.value
        rV = np.sqrt(xV**2+yV**2+zV**2)*u.kpc
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(rV <= rv_max)
        
        # determine the velocity and mass of those particles within the mas radius
        # write your own code below
        vx_new =  self.vx[indexV]
        vy_new =  self.vy[indexV]
        vz_new =  self.vz[indexV]
        m_new =   self.m[indexV]
        
        # compute the center of mass velocity using those particles
        # write your own code below
        vx_COM, vy_COM, vz_COM = self.COMdefine(vx_new, vy_new, vz_new, m_new)
        
        # create an array to store the COM velocity
        # write your own code below
        v_COM = np.array([vx_COM, vy_COM, vz_COM])

        # return the COM vector
        # set the correct units usint astropy
        # round all values                                                                                        
        v_COM = np.round(v_COM,2)
        v_COM =  v_COM*u.km/u.s                 # assinging units
        return v_COM