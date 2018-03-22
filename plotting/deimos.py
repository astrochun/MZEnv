import numpy as np

import matplotlib.patches as mpatches
from matplotlib.path import Path as Path0

from astropy import log

def in_deimos_field(tab0, verts0, silent=False, verbose=True):
    '''
    Determine sources in Keck/DEIMOS field

    Parameters
    ----------
    tab0 : astropy.table.table.Table
      Astropy Table of HSC-SSP NB excess emitter catalog.

    verts0 : list
      List of matplotlib.patches.Rectangle vertices

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
    in_fields0 : list
      List of numpy arrays

    Notes
    -----
    Created by Chun Ly, 3 March 2018

    Modified by Chun Ly, 4 March 2018
     - Bug fix: Incorrect variable name, patches -> verts0
    '''

    if silent == False: log.info('### Begin in_deimos_field : '+systime())

    ra  = tab0['ra'].data
    dec = tab0['dec'].data
    coords = np.array([ra, dec]).transpose()

    n_ptgs = len(verts0)

    in_field0 = []
    for cc in range(n_ptgs):
        t_path = Path0(verts0[cc])

        cp_res = t_path.contains_points(coords)
        in_field = np.array([xx for xx in range(len(tab0)) if cp_res[xx] == True])
        print cc, len(in_field)
        in_field0.append(in_field)
    #endfor

    if silent == False: log.info('### End in_deimos_field : '+systime())
    return in_field0

#enddef


def plot_deimos_fov(ax, coord, maskno, pa=0.0):
    '''
    Overlay Keck/DEIMOS FoV on image

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
      Matplotlib axes from plt.subplots()

    coord : list of list
      Central RA,Dec coordinate provided in degrees. Each position is provided
      as a list within the main list.
      coord = [ [RA1,Dec1], [RA2,Dec2], ...]

    maskno : list
      List of strings for mask name shorthand to overlay on top of mask

    pa : list
      List of PA in degrees for DEIMOS orientation. Positive is E of North.
      Default: 0.0 deg

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
     ax with Rectangular patch added

     verts0 : list
       List containing vertices for each Rectangle patch

    Notes
    -----
    Created by Chun Ly, 1 March 2018
     - Bug fix: typos, return ax
    Modified by Chun Ly, 3 March 2018
     - Loop over each position to overlay Rectangle patch
     - Sign handling and matrix rotation fix so positive PA is indeed E of N
     - Update documentation
     - Get vertices of Rectangle patch and return list of vertices
    Modified by Chun Ly, 5 March 2018
     - Add maskno input for ax.text annotation
    Modified by Chun Ly, 9 March 2018
     - Change color/style scheme for mpatches.Rectangle
    '''

    dx0 =  5.0 / 60.0 # in deg
    dy0 = 16.7 / 60.0 # in deg

    xbox0 = np.array([ 1, 1,-1,-1])*dx0/2.0
    ybox0 = np.array([-1, 1, 1,-1])*dy0/2.0

    n_ptgs = len(coord)

    verts0 = []

    for cc in range(n_ptgs):
        rad = -1*pa[cc]*np.pi/180.0 #Sign flip

        xprime =  xbox0*np.cos(rad) - ybox0*np.sin(rad) #Sign fix
        yprime =  xbox0*np.sin(rad) + ybox0*np.cos(rad) #Sign fix

        xbox = coord[cc][0] + xprime
        ybox = coord[cc][1] + yprime

        xy_low = [xbox[3], ybox[3]]

        t_rot0 = -90+pa[cc] if pa[cc] > 0 else 90+pa[cc]
        ax.text(coord[cc][0], coord[cc][1], maskno[cc], ha='center', va='center',
                rotation=t_rot0)

        rect = mpatches.Rectangle(xy_low, dx0, dy0, angle=-1*pa[cc],
                                  linewidth=1.5, ec="k", color="none") # Sign flip

        # Must be done before ax.add_patch call
        verts0.append(rect.get_verts())
        ax.add_patch(rect)
    #endfor

    return ax, verts0
#enddef
