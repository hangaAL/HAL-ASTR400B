# imports
from photutils.isophote import EllipseGeometry
from photutils.isophote import Ellipse
import matplotlib.pyplot as plt

# functions to fit ellipses using photutils.isophote

# based off code from ellipse_fitting Jupyter notebook by Himansh

def fit_ellipse_isophote(hist, sma_base, eps_base, pa_base):
    '''
    Fits elliptical isphotes to the given 2D  galaxy particle 
    density histogram using `photutils.isophote`.

    Parameters
    ----------
        hist: `numpy.histogram2d`
            2D particle density histogram
        sma_base: `float`
            starting point for the semimajor
            axis fit
        eps_base: `float`
            starting point for the 
            ellipticity fit
        pa_base: `float`
            starting point for the position
            angle fit

    Returns
    -------
        ell_list: `photutils.isophote.IsophoteList`
            the elliptical isophote fits for the
            given galaxy.
    '''

    # unpack particle density histogram & x & y bin edges
    dens, xedges, yedges = hist[0], hist[1], hist[2]
    # get x & y bin sizes
    xbinsize = xedges[1] - xedges[0]
    ybinsize = yedges[1] - yedges[0]

    # get pixel coordinates of center
    center_x_pix = -xedges[0]/xbinsize
    center_y_pix = -yedges[0]/ybinsize

    # define base ellipse geometry for fitting
    # sma is semimajor axis, eps is ellipticity, pa is position angle
    ell_geom = EllipseGeometry(x0 = center_x_pix, y0 = center_y_pix,
                               sma=sma_base/xbinsize, eps=eps_base, 
                               pa=pa_base)
    # initialize the ellipse object for fitting
    # transpose the histogram to undo the transposition done by Ellipse
    ellipse_base = Ellipse(dens.T, ell_geom)
    # fit the ellipses
    ell_list = ellipse_base.fit_image()

    return ell_list


# also based off code from ellipse_fitting Jupyter notebook by Himansh

def get_closest_ellipses(sma_arr, ell_list, binsize, ax, **kwargs):
    '''
    Finds the elliptical isophotes with semimajor axes closest 
    to those in the given array. 

    Parameters
    ----------
        sma_arr: `numpy.ndarray-like`
            the array of target semimajor
            axes
        ell_list: `photutils.isophote.IsophoteList`
            the elliptical isophote fits for the
            given galaxy
        binsize: `float`
            the binsize of the projection, 
            used for converting from pixels
            back to kpc
        ax: `matplotlib.axes.Axes`
            the axis (subplot) to place the plotted
            ellipses on 
        **kwargs:
            arguments for plotting the ellipses
            (color, linestyle, etc.) 

    Returns
    -------
        closest_ells: `list`
            list of the isophotes with semimajor
            axes closest to the targets
    '''
    # convert semimajor axes to pixels
    sma_pix = sma_arr/binsize

    # using a list here because isophotes don't lend themselves
    # well to being put in an array
    # get isophotes with semi-major axes closest to the given ones
    closest_ells = [ell_list.get_closest(sma) for sma in sma_pix]

    # get x, y points associated with these ellipses and plot them
    for ellipse in closest_ells:
        x_ell, y_ell = ellipse.sampled_coordinates()
        ax.plot(x_ell, y_ell, **kwargs)
    
    return closest_ells
