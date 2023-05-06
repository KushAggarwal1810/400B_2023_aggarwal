"""

Astr 400 B
Final Project
Kush Aggarwal

"""

"""
The aim of this code to read in the data files both before merger and after merger
and plot the velocity dispersion and bh mass graohs. 

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


"""
Reading data and plotting
"""

file = open('sigma_MW.txt','r' )                
data = np.genfromtxt('sigma_MW.txt',dtype=None,names=True)  # storing all labels to data

t_MW = data['t']
sigx_MW = data['sigx']
sigy_MW= data['sigy']
sigz_MW= data['sigz']
mbx_MW= data['mbx']
mby_MW= data['mby']
mbz_MW= data['mbz']


file = open('sigma_M31.txt','r' )                
data = np.genfromtxt('sigma_M31.txt',dtype=None,names=True)  # storing all labels to data

t_M31 = data['t']
sigx_M31 = data['sigx']
sigy_M31= data['sigy']
sigz_M31= data['sigz']
mbx_M31= data['mbx']
mby_M31= data['mby']
mbz_M31= data['mbz']

file = open('sigma_Mer.txt','r' )                
data = np.genfromtxt('sigma_Mer.txt',dtype=None,names=True)  # storing all labels to data

t_Mer = data['t']
sigx_Mer = data['sigx']
sigy_Mer= data['sigy']
sigz_Mer= data['sigz']
mbx_Mer= data['mbx']
mby_Mer= data['mby']
mbz_Mer= data['mbz']


t1= np.concatenate((t_MW,t_Mer))


"MW net data before and after merger"

sig_xMW= np.concatenate((sigx_MW,sigx_Mer))
sig_yMW= np.concatenate((sigy_MW,sigy_Mer))
sig_zMW= np.concatenate((sigz_MW,sigz_Mer))
mb_xMW= np.concatenate((mbx_MW,mbx_Mer))
mb_yMW= np.concatenate((mby_MW,mby_Mer))
mb_zMW= np.concatenate((mbz_MW,mbz_Mer))

"M31 net data before and after merger"

sig_xM31= np.concatenate((sigx_M31,sigx_Mer))
sig_yM31= np.concatenate((sigy_M31,sigy_Mer))
sig_zM31= np.concatenate((sigz_M31,sigz_Mer))
mb_xM31= np.concatenate((mbx_M31,mbx_Mer))
mb_yM31= np.concatenate((mby_M31,mby_Mer))
mb_zM31= np.concatenate((mbz_M31,mbz_Mer))




# plt.plot(t1,sig_xMW)
# plt.plot(t1,sig_xM31)
# # # plt.plot(t1,sig_yMW)
# # # plt.plot(t1,sig_zMW)
# plt.show()

# plt.plot(t1,mb_xMW)
# plt.plot(t1,mb_xM31)

i= np.where(mb_xMW ==min(mb_xMW))
j= np.where(mb_xM31 ==min(mb_xM31))

print(len(mbx_M31))
print(len(mbx_MW))

print(mbx_Mer[0])
# print(j)

# print(mb_xMW[88])
# print(mb_xMW[89])

# print(mb_xM31[88])
# print(mb_xM31[89])


fig = plt.figure()
gs = fig.add_gridspec(3, hspace=0)
axs = gs.subplots(sharex=True, sharey=True)
fig.suptitle('Velocity dispersion for MW')
axs[0].plot(t1, sig_xMW, 'red',label=r'$\sigma_x$' )
axs[0].legend(loc='upper left', fontsize='12')
axs[1].plot(t1, sig_yMW, 'blue',label= r'$\sigma_y$')
axs[1].legend(loc='upper left', fontsize='12')
axs[2].plot(t1, sig_zMW, 'green',label= r'$\sigma_z$')
axs[2].legend(loc='upper left', fontsize='12')

# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
fig.supxlabel('Time (Gyr)')
fig.supylabel('Velocity dispersion (km/s)')
plt.show()

fig = plt.figure()
gs = fig.add_gridspec(3, hspace=0)
axs = gs.subplots(sharex=True, sharey=True)
fig.suptitle('Velocity dispersion for M31')
axs[0].plot(t1, sig_xM31, 'red',label=r'$\sigma_x$' )
axs[0].legend(loc='upper left', fontsize='12')
axs[1].plot(t1, sig_yM31, 'blue',label= r'$\sigma_y$')
axs[1].legend(loc='upper left', fontsize='12')
axs[2].plot(t1, sig_zM31, 'green',label= r'$\sigma_z$')
axs[2].legend(loc='upper left', fontsize='12')

# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
fig.supxlabel('Time (Gyr)')
fig.supylabel('Velocity dispersion ( km/s)')
plt.show()


fig = plt.figure()
gs = fig.add_gridspec(3, hspace=0)
axs = gs.subplots(sharex=True, sharey=True)
fig.suptitle('Black Hole Mass for MW')
axs[0].plot(t1, mb_xMW, 'red',label=r'$BM_x$' )
axs[0].legend(loc='upper left', fontsize='12')
axs[1].plot(t1, mb_yMW, 'blue',label= r'$BM_y$')
axs[1].legend(loc='upper left', fontsize='12')
axs[2].plot(t1, mb_zMW, 'green',label= r'$BM_z$')
axs[2].legend(loc='upper left', fontsize='12')

# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
fig.supxlabel('Time (Gyr)')
fig.supylabel('Black Hole Mass ($10^{7} M_{\odot}$)')
plt.show()

fig = plt.figure()
gs = fig.add_gridspec(3, hspace=0)
axs = gs.subplots(sharex=True, sharey=True)
fig.suptitle('Black Hole Mass for M31')
axs[0].plot(t1, mb_xM31, 'red',label=r'$BM_x$' )
axs[0].legend(loc='upper left', fontsize='12')
axs[1].plot(t1, mb_yM31, 'blue',label= r'$BM_y$')
axs[1].legend(loc='upper left', fontsize='12')
axs[2].plot(t1, mb_zM31, 'green',label= r'$BM_z$')
axs[2].legend(loc='upper left', fontsize='12')

# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
fig.supxlabel('Time (Gyr)')
fig.supylabel('Black Hole Mass ($10^{7} M_{\odot}$)')
plt.show()









