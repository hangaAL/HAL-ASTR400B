import numpy as np
import astropy.units as u

from ReadFile import Read
from CenterOfMass import CenterOfMass
from astropy.constants import G

class MassProfile:

    def __init__(self, galaxy, snap, extpath=None):
        '''
        Class to determine the mass profile and rotation curves 
        for a galaxy at a specific snapshot.

        -----------------
        Parameters:
            galaxy (str): 
                the name of the galaxy to use
            snap (int):
                the snap number of the galaxy
            extpath (str or None; default=None):
                path to the galaxy file if it's
                not in the current directory
        '''
        
        # add the filenumber to the value "000" in the filename
        ilbl = '000' + str(snap)
        # remove all but the last 3 digits for consistency
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy) + ilbl + '.txt'
        if extpath != None:
            self.filename = extpath + '/' + self.filename

        __, __, data = Read(self.filename) # get galaxy data

        # get positions, masses, and particle types
        self.x = data['x']*u.kpc
        self.y = data['y']*u.kpc
        self.z = data['z']*u.kpc
        self.mass = data['m']
        self.particles = data['type']

        # get G in the units we need
        # also this is the only way it works for some reason
        self.G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

        # get galaxy name
        self.gname = galaxy

    def MassEnclosed(self, ptype, r):
        '''
        Calculates the enclosed mass at the specified radii
        for a given particle type.

        ---------------
        Parameters:
            ptype (int; 1, 2 or 3):
                the type of particle to look at
            r (numpy array of floats):
                the radii to calculate the 
                enclosed mass at
        
        ---------------
        Returns:
            m_r (numpy array of astropy quantities):
                the mass profile for that particle 
                in M_sun

        '''

        # get COM position
        COM = CenterOfMass(self.filename, 2)
        COM_p = COM.COM_P()
        
        # get indices corresponding to the particle type
        p_idx = np.where(self.particles == ptype)
        # get positions of just those particles
        x_new = self.x - COM_p[0]
        y_new = self.y - COM_p[1]
        z_new = self.z - COM_p[2]
        # calculate radii of the particles
        r_all = np.sqrt(x_new**2 + y_new**2 + z_new**2)

        # make sure r has the right units
        r = r*u.kpc
        # make an array to save all the masses to
        nr = len(r)
        m_r = np.zeros(nr)

        # loop through all the radii
        for i in range(nr):

            # get indices of particles within the current radius
            r_idx = np.where((r_all <= r[i])&(self.particles == ptype))
            # sum their masses
            mtot = np.sum(self.mass[r_idx])
            # save it to the array
            m_r[i] = mtot

        m_r = m_r*1e10*u.Msun # give it the right units

        return m_r
    
    def MassEnclosedTotal(self, r):
        '''
        Calculates the enclosed mass at the specified radii
        for all particle types.

        ---------------
        Parameters:
            r (numpy array of floats):
                the radii to calculate the 
                enclosed mass at
        
        ---------------
        Returns:
            mtot (numpy array of astropy quantities):
                the total mass profile in M_sun

        '''
        
        # get mass profiles for each particle type
        m_halo = self.MassEnclosed(1, r)
        m_disk = self.MassEnclosed(2, r)
        m_bulge = self.MassEnclosed(3, r)

        # sum them to get total mass
        mtot = m_halo + m_disk + m_bulge

        return mtot
    
    def HernquistMass(self, r, a, Mhalo):
        '''
        Calculates the Hernquist mass profile for a
        given halo mass and a.

        -----------
        Parameters:
            r (numpy array of floats):
                the radii to calculate the 
                enclosed mass at
            a (astropy quantity):
                scale radius in kpc
            Mhalo (astropy quantity):
                total dark matter mass in Msun

        ------------
        Returns:
            M_hern (numpy array of astropy quantities):
                Hernquist mass profile
        '''


        r = r*u.kpc
        # numerator
        num = Mhalo*r**2
        # denominator
        denom = (a + r)**2

        # calculate Hernquist profile
        M_hern = (num/denom)

        return M_hern
    
    def CircularVelocity(self, ptype, r):
        '''
        Calculates the rotation curve for the specified radii
        for a given particle type.

        ---------------
        Parameters:
            ptype (int; 1, 2 or 3):
                the type of particle to look at
            r (numpy array of floats):
                the radii to calculate the 
                enclosed mass at
        
        ---------------
        Returns:
            vcirc_r (numpy array of astropy quantities):
                the rotation curve for that particle 
                in M_sun

        '''
        # make array to store velocities
        nr = len(r)
        vcirc_r = np.zeros(nr)*u.km/u.s

        # get mass profile of particle
        mr = self.MassEnclosed(ptype, r)

        # loop through radii
        # and calculate circular velocity at each one
        for i in range(nr):
            
            rad = r[i]*u.kpc
            M = mr[i]
            vc = np.sqrt(self.G*M/rad)
            vcirc_r[i] = np.around(vc, 2)

        return vcirc_r
    
    def CircularVelocityTotal(self, r):
        '''
        Calculates the rotation curve for the specified radii
        for all particle type.

        ---------------
        Parameters:
            r (numpy array of floats):
                the radii to calculate the 
                enclosed mass at
        
        ---------------
        Returns:
            vcirctot_r (numpy array of astropy quantities):
                the total rotation curve in M_sun

        '''
        
        # make array to store velocities
        nr = len(r)
        vcirctot_r = np.zeros(nr)*u.km/u.s

        # get total mass profile 
        mtot_r = self.MassEnclosedTotal(r)

        # loop through radii
        # and calculate circular velocity at each one
        for i in range(nr):
            
            rad = r[i]*u.kpc
            M = mtot_r[i]
            vc = np.sqrt(self.G*M/rad)
            vcirctot_r[i] = np.around(vc, 2)

        return vcirctot_r
    
    def HernquistVCirc(self, r, a, Mhalo):
        '''
        Calculates the Hernquist rotation curve for a
        given halo mass and a.

        -----------
        Parameters:
            r (numpy array of floats):
                the radii to calculate the 
                enclosed mass at
            a (astropy quantity):
                scale radius in kpc
            Mhalo (astropy quantity):
                total dark matter mass in Msun

        ------------
        Returns:
            vc_hern (numpy array of astropy quantities):
                Hernquist rotation curve
        '''
        
        # calculate Hernquist mass profile
        M_hern =  self.HernquistMass(r, a, Mhalo)
        r = r*u.kpc
        # calculate rotation curve for mass profile
        vc_hern = np.around(np.sqrt(self.G*M_hern/r), 2)

        return vc_hern