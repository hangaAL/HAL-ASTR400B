# imports
import numpy as np
import astropy.units as u
from ReadFile import Read

def ComponentMass(fname, component):
    '''
    Calculates the total mass of a component of a galaxy for
    the given file.

    Parameters:
        fname (string): the path to the file containing the 
            galaxy data.
        component (string): the galaxy component whose mass
            is to be calculated.
    
    Returns:
        total_mass (astropy quantity): the total mass of the 
        galaxy component, in units of 10^12 Msun.
    '''

    ptype_dict = {"Halo":1., "Disk":2., "Bulge":3.} # assigning the correct type number to each type
    ptype_num = ptype_dict[component] # getting the type number of the input component

    __, __, data = Read(fname) # read in the particle data; time and total particle number not needed
    data_filtered = data[data['type']==ptype_num] # filtering the data for particles of the right type

    total_mass = np.around(np.sum(data_filtered['m'])/100, 3)*u.M_sun #originally in 10^10 Msun, 
    #so divide by 100 to get it in 10^12 Msun units. also rounded to 3 decimal places

    return total_mass