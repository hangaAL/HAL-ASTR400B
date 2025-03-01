

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




def OrbitCOM(galaxy, start, end, n):
    """
    Loops over all the desired snapshots to compute 
    the COM pos and vel as a function of time.

    Parameters:
        galaxy (str):
            the name of the galaxy
        start (int):
            the number of the first snapshot to be read in
        end (int):
            the number of the last snapshot to be read in
        n (int):
            the interval over which the COM will be returned
          
    Returns:
        None 
    """
    
    # compose the filename for output
    fileout = 'Orbit_' + galaxy + '.txt'
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1
    if galaxy == 'M33': # M33 will be severely tidally disrupted
        volDec = 4
    else:
        volDec = 2
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start, end+n, n)
    if snap_ids.size == 0:
        sys.exit("Please enter a valid snapshot ID range.")
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    n_snap = snap_ids.size # number of snap IDs
    orbit = np.zeros((n_snap, 7))
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over files
        
        # compose the data filename (be careful about the folder)
        ilbl = '000' + str(snap_id)
        filein = galaxy + '/' + galaxy + '_' + ilbl[-3:] +'.txt' # make sure snap ID has 3 digits
        # Initialize an instance of CenterOfMass class, using disk particles
        gal_COM = CenterOfMass(filein, 2)
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        gCOM_p = gal_COM.COM_P(volDec)
        gCOM_v = gal_COM.COM_V(gCOM_p[0], gCOM_p[1], gCOM_p[2])
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        orbit[i] = (gal_COM.time.value/1000, *tuple(gCOM_p.value), 
                    *tuple(gCOM_v.value))
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
    
# SEE HOMEWORK 6 JUPYTER NOTEBOOK FOR ANSWERS