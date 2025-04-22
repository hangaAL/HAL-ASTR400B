# Density contour plotting code from Lab 7

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

def density_contour(xdata, ydata, nbins_x, nbins_y, conf_lvls, ax=None, **contour_kwargs):
    """ Create a density contour plot.
        NOTE: Modified from Lab 7 code to allow for any given confidence
        levels as input.
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    conf_lvls : list or np.ndarray
        Confidence levels to use for contours, in the order 
        in which they are to be stored
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
        
    Example Usage
    -------------
     density_contour(x pos, y pos, contour res, contour res, axis, colors for contours)
     e.g.:
     density_contour(xD, yD, 80, 80, ax=ax, 
         colors=['red','orange', 'yellow', 'orange', 'yellow'])

    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), density=True)
    # NOTE : if you are using the latest version of python, in the above: 
    # instead of normed=True, use density=True
    
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))
    
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T
    fmt = {}
    
    ### Adjust Here #### 
    
    # Contour Levels Definitions

    # to store contour levels
    levels = []

    # fuck, shit = plt.subplots()

    # fu = np.linspace(0,1,100)
    # egg = np.zeros((3,100))
    # for i in range(3):
    #     for j in range(100):
    #         egg[i,j] = find_confidence_interval(fu[j], pdf, conf_lvls[i])
    
    # shit.plot(fu, egg[0])
    # shit.plot(fu, egg[1])
    # shit.plot(fu, egg[0])
    
    # get contour levels
    for conf_lvl in conf_lvls:
        # calculate level
        level = so.brentq(find_confidence_interval, 0., 1., args=(pdf, conf_lvl))
        # and store it
        levels.append(level)

    # contour level labels
    strs = [str(conf_lvl) for conf_lvl in conf_lvls]

    # one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    # two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    # three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    
    # You might need to add a few levels

    # Array of Contour levels. Adjust according to the above
    # levels = [one_sigma, two_sigma, three_sigma][::-1]
    
    # contour level labels  Adjust accoding to the above.
    # strs = ['0.68', '0.95', '0.99'][::-1]

    
    ###### 
    
    if ax == None:
        contour = plt.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        plt.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)

    else:
        contour = ax.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        ax.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)
    
    return contour