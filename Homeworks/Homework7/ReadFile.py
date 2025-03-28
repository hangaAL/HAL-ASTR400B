# Contains the Read() function for reading in
# and outputting particle properties from a galaxy
# simulation at a certain point in time

# Imports
import numpy as np        # for arrays
import astropy.units as u # for units

def Read(filename):
    '''
    Reads in a text file containing particle data for a simulated
    galaxy at a specific time in the simulation, and outputs the 
    particle data with identifying information.

    Parameters:
        filename (str): the path to the file to be read.

    Returns:
        time (astropy quantity): the time in the 
            simulation this particle data corresponds
            to in Myr
        n_tot (int): the total number of particles in
            the simualtion.
        data (numpy ndarray): all of the particle data. 
            Columns are: type, m, x, y, z, vx, vy, vz
    '''

    file = open(filename, "r") # opening the file

    # reading in the time and total number of particles
    label, value = file.readline().split() # getting first line
    time = float(value)*u.Myr # saving time in Myr
    label, value = file.readline().split() # getting second line
    n_tot = int(value) # saving total number of particles

    file.close() # closing the file so it doesn't hog memory

    data = np.genfromtxt(filename, dtype=None, 
                         names=True, skip_header=3) # saving particle data as an ndarray
    
    return time, n_tot, data