import numpy as np

import matplotlib.patches as mpatches

from astropy import log

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
