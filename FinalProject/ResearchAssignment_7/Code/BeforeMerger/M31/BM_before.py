"""

Astr 400 B
Final Project
Kush Aggarwal

"""

"""
The aim of this code to find the black hole mass in the x,y,z directions using the
M-sigma relation

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


def M_sigma(a):
    
    """
    a: array of velocity dispersions in different directions
    """
    # M_sun= 2*10**30
    A= 3.09
    alpha=4.4
    
    M= 10*A*(a/200)**alpha             # ( Kastytis Zubovas and Andrew R King, 2019)
    
    return M                 # M is in  10**7 solar mass units
    
