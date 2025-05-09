# make function to calculate COM of two COMs
def COM_of_COMs(COMvec1, COMvec2, M1, M2):
    ''' 
    Calculates the center of mass position or 
    velocity for two given centers of mass and
    the associated total masses.

    Parameters
    ----------
        COMvec1, COMvec2: `numpy array-like`
            the vectors associated with each
            center of mass respectively
        M1, M2: `astropy quantity`
            the masses associated with each
            center of mass (Msun)
    
    Returns
    -------
        COMvec_tot: `numpy array-like`
            the vector for the total
            center of mass
    '''
    # calculate total center of mass vector
    denom = M1+M2
    num = COMvec1*M1 + COMvec2*M2
    COMvec_tot = num / denom

    return COMvec_tot