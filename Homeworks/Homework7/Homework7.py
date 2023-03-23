

"""
Here we create a series of functions that will determine the acceleration M33
feels from M31 and integrate its current position and velocity forwards in time.

Kush Aggarwal

Astr 400B

Homework 7
"""


# import necessary modules
import numpy as np
import astropy.units as u
import astropy.table as tbl
from numpy import linalg as LA
from ReadFile import Read
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
# from CenterOfMass import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass
from CenterOfMass2 import CenterOfMass   # usign updated center of mass function.
# # M33AnalyticOrbit




class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): # **** add inputs
    
        """ Here we are definig the glbal constants
          inputs : filename 
          Nmae of the file that is the output
        """
                
        ### get the gravitational constant (the value is 4.498502151575286e-06)     
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
           
        ### **** store the output file name
        self.filename = "%s_" % filename + '.txt'
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        M33_COM = CenterOfMass("M33_000.txt", 2)

        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        MW_COM_p_M33 = M33_COM.COM_P(0.1,4)
        
        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        M33_COM_v = M33_COM.COM_V(MW_COM_p_M33[0], MW_COM_p_M33[1], MW_COM_p_M33[2])
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        M31_COM = CenterOfMass("M31_000.txt", 2)
        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        MW_COM_p_M31 = M31_COM.COM_P(0.1,2)
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        M31_COM_v = M31_COM.COM_V(MW_COM_p_M31[0], MW_COM_p_M31[1], MW_COM_p_M31[2])
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0 = (MW_COM_p_M33-MW_COM_p_M31).value
        self.v0 = (M33_COM_v - M31_COM_v).value
        ### get the mass of each component in M31 
        ### disk
        # **** self.rdisk = scale length (no units)
        self.rdisk= 5
        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk= ComponentMass("M31_000.txt", 2)*1e12
        ### bulge
        # **** self.rbulge = set scale length (no units)
        self.rbulge=1
        
        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge= ComponentMass("M31_000.txt", 3)*1e12
        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 62
        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo= ComponentMass("M31_000.txt", 1)*1e12
    
    
    def HernquistAccel_buldge(self,M, r_a,r): # it is easiest if you take as an input the position VECTOR 
        """ 
        
        This function computes the graviatational acceleration due to buldge particles of M31
        
        Inputs :
            
          M : Mass of bulde particles of M31
          r_a = scaling length from homework 5
          r : position vector
         
        Outputs :
            returns the acceleration induced by buldge particles of M31
        
        """
        ### **** Store the magnitude of the position vector
        rmag = LA.norm(r)
         
        ### *** Store the Acceleration  # using hernquist profile
        Hern =  -(self.G*M*r)/(rmag*(r_a+rmag)**2) #follow the formula in the HW instructions
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern
    
    def HernquistAccel_halo(self,M, r_a,r): # it is easiest if you take as an input the position VECTOR 
        """ 
        
        This function computes the graviatational acceleration due to HALO PARTICLES of M31
        
        Inputs :
            
          M : Mass of halo particles of M31
          r_a = scaling length from homework 5
          r : position vector
         
        Outputs :
            returns the acceleration induced by halo particles of M31 """
    
        ### **** Store the magnitude of the position vector
        rmag = LA.norm(r)
        
        ### *** Store the Acceleration    # using hernquist profile
        Hern =  -(self.G*M*r)/(rmag*(r_a+rmag)**2) #follow the formula in the HW instructions
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern
    
    def MiyamotoNagaiAccel(self,M,r_d,r):# it is easiest if you take as an input a position VECTOR  r 
        """ 
        
        This function computes the graviatational acceleration due to disk particles of M31
        
        Inputs :
            
          M : Mass of disk particles of M31
          r_a = scaling length from homework 5
          r : position vector
         
        Outputs :
            returns the acceleration induced by disk particles of M31 """
        z_d = r_d/5
        R= np.sqrt(r[0]**2 + r[1]**2)
        B= r_d+ np.sqrt(r[2]**2+ z_d**2)
        
        
        ### Acceleration **** follow the formula in the HW instructions
        m = np.array([1,1,B/(np.sqrt(r[2]**2+z_d**2))])
        
        a = -(self.G*M*r*m)/(R**2 +B**2)**(1.5)   # Miyamoto-Nagai 1975 profile
        
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        
        
       
        return a
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self,r): # input should include the position vector, r
        """ **** ADD COMMENTS """
        
        ### Call the previous functions for the halo, bulge and disk
        x=self.HernquistAccel_buldge(self.Mbulge,self.rbulge, r)
        y= self.HernquistAccel_halo(self.Mhalo,self.rhalo,r)
        z= self.MiyamotoNagaiAccel(self.Mdisk,self.rdisk,r)
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        a_total= x+y+z   
            # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return a_total
    
    
    
    def LeapFrog(self,dt,r,v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """ **** ADD COMMENTS """
        
        # predict the position at the next half timestep
        rhalf = r+v*dt/2
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        a= self.M31Accel(rhalf).value    # finding acceleration
        # print(a)
        vnew = v+ a*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf + vnew*dt/2
        
        return  rnew, vnew # **** return the new position and velcoity vectors
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """ we will loop over the LeapFrog integrator to solve the equations of motion and compute
the future orbit of M33 for 10 Gyr into the future.

        Inputs :
        
  t0: intial time
    dt : time step
tmax : final time            



        """

        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros((int(tmax/dt)+2 ,7) )
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        pos=np.array([])
        vel= np.array([])
        pos= np.append(pos,self.r0)  # initialising the position and velocity
        vel= np.append(vel,self.v0)
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t <= tmax):  # as long as t has not exceeded the maximal time 
            # print(pos,vel)
            
            # **** advance the time by one timestep, dt
            t= t+dt
            # **** store the new time in the first column of the ith row
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            
            a,b = self.LeapFrog(dt, pos,vel)
             
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
           
            
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            orbit[i] = t , *tuple(a), *tuple(b)
            
            pos=np.array([])     
            vel= np.array([])
            pos= np.append(pos,a)
            vel=np.append(vel,b)
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i= i+1
        
        
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                    header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                    .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function

if __name__ == '__main__' : 
     
      M= M33AnalyticOrbit('orbit2')                 # creating an instance of class
      b= M.OrbitIntegration(0, 0.1,10)              # calling the orbit integrator with inital time 0, time step 0.1 and
                                                    # final time 10 Gyr


"""
The above code should be run first, it will generate a orbit file, now comment the code above this line 
and then run the below code
"""    
file = open('orbit2_.txt','r' )       # openign created file                   
data = np.genfromtxt('orbit2_.txt',dtype=None,names=True)  # storing all labels to data

t_MW = data['t']
x_MW = data['x']
y_MW= data['y']
z_MW= data['z']
vx_MW= data['vx']
vy_MW= data['vy']
vz_MW= data['vz']

rel_sep= np.array([])
rel_vel= np.array([])

for i in range(len(t_MW)):
    a = np.array([x_MW[i],y_MW[i],z_MW[i]])
    b=  np.array([vx_MW[i],vy_MW[i],vz_MW[i]])
    v=LA.norm(a)
    w= LA.norm(b)
    rel_sep= np.append(rel_sep,v)
    rel_vel= np.append(rel_vel,w)


"""
Copying the code from Homework 6, so as to overplot the soltuoj from this code.
"""
    
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

def mag(a,b):
    diff = np.array([])
    for i in range(len(a)):
        t= abs(a[i]-b[i])
        diff = np.append(diff,t)
    v=LA.norm(diff)
    return v   

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
    
"""
The code betweent this comment and last comment has been copied from hw 6
"""

plt.plot(t_MW, rel_sep,color='green',label=r'M33-M31(homework7)')
plt.plot(t_MW, rel_sep3,color='red',label=r'M33-M31 (homework6)')
plt.xlabel('Time(Gyr)')
plt.ylabel('Separation(kpc)')
plt.title('Figure for Relative separation between M33 and M31')
plt.legend()
plt.show()

plt.plot(t_MW, rel_vel,color='green',label=r'M33-M31(homework7)')
plt.plot(t_MW, rel_sep4,color='blue',label=r'M33-M33 (homework6)')
plt.xlabel('Time(Gyr)')
plt.ylabel('velocity(kpc)')
plt.title('Figure for Relative velocity between M33 and M33')
plt.legend()
plt.show()