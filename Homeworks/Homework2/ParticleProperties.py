import numpy as np
import astropy.units as u
from ReadFile import Read

def ParticleInfo(filename, ptype, pnum):
    '''
    Finds the specified particle in the given data file and
    returns the magnitude of its distance from the center of
    mass of the Milky Way and its velocity, along with its 
    mass.

    Parameters:
        filename (string): the path to the file to be read.
        ptype (string): the type of the particle. Can be
            Halo, Disk, or Bulge.
        pnum (int): the order of the particle within its
            type (i.e. pnum=5 for the 5th Halo particle)

    Returns:
        dist (astropy quantity): the magnitude of the 
            distance of the particle from the Milky 
            Way center of mass (kpc).
        vel (astropy quantity): the magnitude of the 
            velocity of the particle (km/s).
        mass (astropy quantity): the mass of the 
            particle (Msun).
    '''

    ptype_dict = {"Halo":1., "Disk":2., "Bulge":3.} # assigning the correct type number to each type
    ptype_num = ptype_dict[ptype] # getting the type number of the input particle

    time, n_tot, data = Read(filename) # getting all data from the file

    data_filtered = data[data['type']==ptype_num] # filtering the data for particles of the right type
    particle = data_filtered[pnum-1] # getting the right particle

    # getting the x, y, & z components of the distance
    x = particle['x']*u.kpc 
    y = particle['y']*u.kpc
    z = particle['z']*u.kpc
    dist = np.around(np.sqrt(x**2 + y**2 + z**2), 3) # the magnitude of the distance, rounded to 3 decimal places

    # getting the x, y, & z components of the velocity
    vx = particle['vx']*u.km/u.s
    vy = particle['vy']*u.km/u.s
    vz = particle['vz']*u.km/u.s
    vel = np.around(np.sqrt(vx**2 + vy**2 + vz**2), 3) # the magnitude of the velocity, rounded to 3 decimal places

    mass = particle['m']*1e10*u.M_sun # the mass of the particle 

    return dist, vel, mass