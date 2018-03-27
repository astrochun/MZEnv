import numpy as np

from chun_codes import systime
import matplotlib.patches as mpatches

from astropy import log

def in_hecto_field(tab0, fld_coord, silent=False, verbose=True):
    '''
    Determine sources in MMT/Hectospec field

    Parameters
    ----------
    tab0 : astropy.table.table.Table
      Astropy Table of HSC-SSP NB excess emitter catalog.

    fld_coord : list
      List of RA,Dec coordinate set

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
    Created by Chun Ly, 16 March 2018
    '''

    if silent == False: log.info('### Begin in_hecto_field : '+systime())

    ra  = tab0['ra'].data
    dec = tab0['dec'].data

    n_ptgs = len(fld_coord)

    in_field0 = []
    for cc in range(n_ptgs):
        ra_diff0  = ra  - fld_coord[cc][0]
        dec_diff0 = dec - fld_coord[cc][1]
        diff0 = np.sqrt(ra_diff0**2 + dec_diff0**2)
        in_field = np.array([xx for xx in range(len(tab0)) if
                             diff0[xx] <= 0.50])
        print cc, len(in_field)
        in_field0.append(in_field)
    #endfor

    if silent == False: log.info('### End in_hecto_field : '+systime())
    return in_field0
#enddef

def plot_hecto_fov(ax, coord, configno):
    '''
    Overlay MMT/Hecto FoV on image

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
      Matplotlib axes from plt.subplots()

    coord : list
      Central RA,Dec coordinate provided in degrees

    Returns
    -------
     ax with Circle patch added

    Notes
    -----
    Created by Chun Ly, 1 March 2018
     - Bug fix for no PA input
     - Aesthetic changes (unfilled circle)
    Modified by Chun Ly, 16 March 2018
     - Loop over each input coordinate set and add Circle patch for each
       Hecto pointing
    Modified by Chun Ly, 27 March 2018
     - Label Hecto pointings with configno
    '''

    radius = 0.5 # in deg

    n_ptgs = len(coord)

    for cc in range(n_ptgs):
        circ = mpatches.Circle(coord[cc], radius=radius, alpha=0.5, ec="b",
                               color="none", linewidth=1.5)
        ax.add_patch(circ)
        ax.text(coord[cc][0], coord[cc][1]-0.475, configno[cc], ha='center',
                va='center', rotation=0.0, color='blue')

    return ax
#enddef
