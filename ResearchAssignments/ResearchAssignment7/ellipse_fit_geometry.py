# imports
import numpy as np
from photutils.isophote import EllipseGeometry
from photutils.isophote import Ellipse

# functions dealing with the geometry of fitted ellipses

# define function to get all the semimajor axes, 
# semiminor axes, position angles, and ellipticities 
# for a projection to make my life easier

def get_geom_properties(ellfits):
    '''
    Gets the geometrical properties of all
    fitted ellipses in a projection.

    Parameters
    ----------
        ellfits: `list`
            the fitted ellipses
    
    Returns
    -------
        smaa_arr: `numpy.ndarray`
            the semimajor axes of
            the fitted ellipses
        smia_arr: `numpy.ndarray`
            the semiminor axes of
            the fitted ellipses
        pa_arr: `numpy.ndarray`
            the position angles of
            the fitted ellipses
        eps_arr: `numpy.ndarray`
            the ellipticities of
            the fitted ellipses
    '''

    n = len(ellfits)
    # storage arrays
    smaa_arr = np.zeros(n) # semimajor axis
    smia_arr = np.zeros(n) # semiminor axis
    pa_arr = np.zeros(n) # position angle
    eps_arr = np.zeros(n) # ellipticity

    # loop through all fitted ellipses
    for i in range(n):

        ellipse = ellfits[i]

        # get and save geometrical properties
        smaa_arr[i] = ellipse.sma
        pa_arr[i] = ellipse.pa
        eps_arr[i] = ellipse.eps
        # from eps = 1 - b/a:
        smia_arr[i] = ellipse.sma*(1-ellipse.eps) 

    return smaa_arr, smia_arr, pa_arr, eps_arr


def get_specific_axis_lengths(geom_params, binsize):
    '''
    Determines which of the ellipses' semimajor
    and semiminor axes is closest to each coordinate
    axis for a set of ellipse geometric parameters.

    Parameters
    ----------
        geom_params: `numpy.ndarray`
            the geometric parameters of the
            fitted ellipses.
        binsize: `float`
            the binsize of the projection, 
            used for converting from pixels
            back to kpc
    
    Outputs
    -------
        spec_axes: `numpy.ndarray`
            the axes sorted into columns
            depending on which coordinate 
            axis they're closest to
            (0 = horizontal, 1 = vertical)
            (ex: xy projection -> 0 = x, 1 = y)
    '''
    # semimajor and minor axes, in kpc
    smaas, smias = geom_params[0]*binsize, geom_params[1]*binsize
    # position angles of smaa from horizontal axis
    pas = geom_params[2]

    # to store the specific axes
    n = len(smaas)
    horiz, vert = np.zeros(n), np.zeros(n)

    # sort them using position angle
    for i in range(n):

        if pas[i] < np.pi/4: # if pa > pi/4, smaa closer to the horizontal
            horiz[i] = smaas[i]
            vert[i] = smias[i]

        else:
            horiz[i] = smias[i]
            vert[i] = smaas[i]

    # condense into one array: 0 = horizontal, 1 = vertical
    spec_axes = np.vstack((horiz, vert))

    return spec_axes


def get_3D_axis_length(axlens, thetas, whichprojs):
    '''
    Calculates an average 3D axis length of an ellipsoid from its 
    projected fits onto 2 planes (e.g. the 3D x-axis
    from its xy and xz projections).

    Parameters
    ----------
        axlens: `numpy.ndarray-like`
            the 2D axis lengths
        thetas: `numpy.ndarray-like`
            the 2D position angles
        whichprojs: `list`
            which 2D axes are being averaged
            (e.g. ['xy', 'yz'])

    Returns
    -------
        len_ax: `float`
            the length of the 3D axis
    '''

    # so I put components in the right part of the vector
    coords_idxs = {'x':0, 'y':1, 'z':2}

    # for storing the vector forms of the 2D axes
    ax_vecs = np.zeros((2,3))

    for i in range(2):
        # get 2D axis vector components
        comp1 = axlens[i]*np.cos(thetas[i])
        comp2 = axlens[i]*np.sin(thetas[i])

        # figure out where to put them using the given
        # projections
        coord1, coord2 = whichprojs[i][0], whichprojs[i][1]
        idx1, idx2 = coords_idxs[coord1], coords_idxs[coord2]

        # make vectors
        ax_vecs[i,idx1] = comp1
        ax_vecs[i,idx2] = comp2
    
    # sum vector components 
    final_ax_vec = np.sum(ax_vecs, axis=0)
    # figure out which axis the two projections have in common using sets
    common_comp = ''.join(set(whichprojs[0]).intersection(whichprojs[1]))
    common_idx = coords_idxs[common_comp]
    # then average that component to get the final 3D axis vector
    final_ax_vec[common_idx] = final_ax_vec[common_idx]/2

    # then get the magnitude of that vector
    len_ax = np.sqrt(np.sum(final_ax_vec**2))

    return len_ax


def triaxiality_param(axes):
    '''
    Calculates the triaxiality parameter of the 
    fitted ellipsoids based on their axis lengths.
    The triaxiality parameter 
    T = (1 - q1^2)/(1 - q2^2) is taken from 
    Quenneville et al. 2022 ApJ 926 30

    Parameters
    ----------
        axes: `np.ndarray`
            the three 3D axes of the 
            fitted ellipsoids.
    
    Returns:
        T: `np.ndarray`
            the triaxiality parameter
        q1: `np.ndarray`
            the ratio of the medium and
            largest axes (B/A)
        q2: `np.ndarray`
            the ratio of the smallest and
            largest axes (C/A)
    '''
    # get largest, medium, and smallest axis lengths
    a = np.max(axes, axis=0)
    b = np.median(axes, axis=0)
    c = np.min(axes, axis=0)

    # calculate axis ratios
    q1, q2 = b/a, c/a

    # get triaxiality parameter
    T = (1 - q1**2)/(1 - q2**2)

    return T, q1, q2