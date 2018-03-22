import numpy as np

from chun_codes import systime
import matplotlib.patches as mpatches
from matplotlib.path import Path as Path0

from astropy import log

def in_bino_field(tab0, verts0, silent=False, verbose=True):
    '''
    Determine sources in MMT/Binospec field

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
    Created by Chun Ly, 22 March 2018
    '''

    if silent == False: log.info('### Begin in_bino_field : '+systime())

    ra  = tab0['ra'].data
    dec = tab0['dec'].data
    coords = np.array([ra, dec]).transpose()

    n_ptgs = len(verts0)

    in_field0 = []
    for cc in range(n_ptgs):
        # Side 1
        t_path1   = Path0(verts0[cc][0])
        cp_res1   = t_path1.contains_points(coords)
        in_field1 = np.array([xx for xx in range(len(tab0)) if
                              cp_res1[xx] == True])

        # Side 2
        t_path2   = Path0(verts0[cc][1])
        cp_res2   = t_path2.contains_points(coords)
        in_field2 = np.array([xx for xx in range(len(tab0)) if
                              cp_res2[xx] == True])

        print cc, len(in_field1), len(in_field2)
        in_field = np.append(in_field1, in_field2)
        in_field.sort()
        in_field0.append(in_field)
    #endfor

    if silent == False: log.info('### End in_bino_field : '+systime())
    return in_field0

#enddef

def plot_bino_fov(ax, coord, maskno, pa=0.0):
    '''
    Overlay MMT/Binospec FoV on image

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
      List of PA in degrees for Binospec orientation. Positive is E of North.
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
    Created by Chun Ly, 21 March 2018
     - Bug fix to handle np.arrays combine and proper lower-left coordinates
     - Simplify handling of ybox1 box corners
     - Bug fixes: Do some np.flip and fix typos: box1 -> boxs1
    '''

    dx0 = 19.2 / 60.0 # for outside; in deg
    dx1 =  3.2 / 60.0 # for inside; in deg
    dy0 = 15.0 / 60.0 # in deg

    xbox0 = np.array([ 1, 1,-1,-1])*dx0/2.0
    xbox1 = np.array([ 1, 1,-1,-1])*dx1/2.0
    ybox0 = np.array([-1, 1, 1,-1])*dy0/2.0

    n_ptgs = len(coord)

    verts0 = []

    for cc in range(n_ptgs):
        rad = -1*pa[cc]*np.pi/180.0 #Sign flip

        # outside
        xprime0 = xbox0*np.cos(rad) - ybox0*np.sin(rad) #Sign fix
        yprime0 = xbox0*np.sin(rad) + ybox0*np.cos(rad) #Sign fix

        # inside
        xprime1 = xbox1*np.cos(rad) - ybox0*np.sin(rad) #Sign fix
        yprime1 = xbox1*np.sin(rad) + ybox0*np.cos(rad) #Sign fix

        # side 1
        xboxs1 = coord[cc][0] + np.append(xprime0[0:2], np.flip(xprime1[0:2], 0))
        yboxs1 = coord[cc][1] + np.append(yprime0[0:2], np.flip(yprime1[0:2], 0))
        # Always get box corner in same order so last entry is lower left

        # side 2
        xboxs2 = coord[cc][0] + np.append(xprime1[2:], xprime0[2:])
        yboxs2 = coord[cc][1] + np.append(np.flip(yprime1[2:],0), yprime0[2:])

        xy_low1 = [xboxs1[3], yboxs1[3]]
        xy_low2 = [xboxs2[3], yboxs2[3]]

        t_rot0 = -90+pa[cc] if pa[cc] > 0 else 90+pa[cc]
        ax.text(coord[cc][0], coord[cc][1], maskno[cc], ha='center',
                va='center', rotation=t_rot0)

        rect1 = mpatches.Rectangle(xy_low1, 8.0/60, dy0, angle=-1*pa[cc],
                                   linewidth=1.5, ec="k", color="none")
        rect2 = mpatches.Rectangle(xy_low2, 8.0/60, dy0, angle=-1*pa[cc],
                                   linewidth=1.5, ec="k", color="none")

        # Must be done before ax.add_patch call
        verts0.append([rect1.get_verts(), rect2.get_verts()])
        ax.add_patch(rect1)
        ax.add_patch(rect2)
    #endfor

    return ax, verts0

#enddef
