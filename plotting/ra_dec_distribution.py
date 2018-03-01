"""
ra_dec_distribution
===================

Plot spatial distribution of HSC-SSP NB excess emitter samples
"""

from chun_codes import systime

from getpass import getuser
from os.path import exists

from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import log

import matplotlib.patches as mpatches # + on 01/03/2018

from .. import read_catalog
from .. import subsample

#bbox_props = dict(boxstyle="square,pad=0.15", fc="white", alpha=0.9, ec="none")

def plot_deimos_fov(ax, coord, pa=0.0):
    '''
    Overlay Keck/DEIMOS FoV on image

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
      Matplotlib axes from plt.subplots()

    coord : list
      Central RA,Dec coordinate provided in degrees

    pa : float
      PA in degrees of DEIMOS orientation. Positive is E of North.
      Default: 0.0 deg

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
     ax with Rectangular patch added

    Notes
    -----
    Created by Chun Ly, 1 March 2018
     - Bug fix: typos, return ax
    '''

    rad = pa*np.pi/180.0

    dx0 =  5.0 / 60.0 # in deg
    dy0 = 16.7 / 60.0 # in deg

    xbox0 = np.array([ 1, 1,-1,-1])*dx0/2.0
    ybox0 = np.array([-1, 1, 1,-1])*dy0/2.0

    xprime =  xbox0*np.cos(rad) + ybox0*np.sin(rad)
    yprime = -xbox0*np.sin(rad) + ybox0*np.cos(rad)
    xbox = coord[0] + xprime
    ybox = coord[1] + yprime

    xy_low = [xbox[3], ybox[3]]
    rect = mpatches.Rectangle(xy_low, dx0, dy0, angle=pa, alpha=0.33, ec="r",
                              color="r")
    ax.add_patch(rect)
    return ax
#enddef

def plot_hecto_fov(ax, coord):
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
    '''

    radius = 0.5 # in deg

    circ = mpatches.Circle(coord, radius=radius, alpha=0.33, ec="b",
                           color="b")
    ax.add_patch(circ)
    return ax
#enddef

def main(field='', dr='pdr1', noOII=False, DEIMOS=False, Hecto=False,
         silent=False, verbose=True):

    '''
    Main function to plot RA and Dec for each sub-sample of NB excess emitters

    Parameters
    ----------
    field : str
      Name of field to read in. Either 'udeep' or 'deep'

    dr : str
      Data release name. Default: 'pdr1'

    noOII : boolean
      Indicate whether to NOT plot [OII] emitters. Default: False

    DEIMOS : boolean
      Indicate whether to plot Keck/DEIMOS FoV on PDF plots

    Hecto : boolean
      Indicate whether to plot MMT/Hecto FoV on PDF plots

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 28 February 2018
     - Import subsample module to get NB excess emitter subsamples
     - Call subsample to get galaxy field, [gal_field]
     - Construct for loop to look over each galaxy field and subsample
     - Define output PDF path
     - Plot ra and dec for each galaxy field and NB subsamples
    Modified by Chun Ly, 1 March 2018
     - Add noOII keyword to prevent [OII] emitters from being plotted
     - log.info output PDF file
     - Add DEIMOS boolean keyword input and plot DEIMOS FoV when set
     - Set axes limit for RA and Dec for each galaxy field
     - ax.legend() visibility improvement, RA/Dec limit changes
     - Bug fix: Placement of n_subsample for ax.legend() ncol determination
     - Add Hecto boolean keyword input and plot Hecto FoV when set
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    if field == '':
        log.warn("### [field] input not specified!!!")
        log.warn("### Either 'udeep' or 'deep'")
        log.warn('### Exiting!!!')
        return

    if getuser() == 'cly':
        dir0 = '/Users/cly/Google Drive/MZEnv/'

    tab0 = read_catalog.main(field=field, dr=dr, silent=silent, verbose=verbose)
    ra0, dec0 = tab0['ra'], tab0['dec']

    sub_dict0, gal_field = subsample.main(tab0=tab0, field=field, dr=dr)

    sample_keys = sub_dict0.keys()

    gal_field0 = list(set(gal_field)) # Get unique list

    RA_lim_dict0 = {'UD-COSMOS':  [149.25,151.00],
                    'UD-SXDS':    [ 33.70, 35.45],
                    'D-COSMOS':   [148.75,151.80],
                    'D-DEEP2_3':  [350.20,353.95],
                    'D-ELAIS_N1': [240.10,245.40]}

    DE_lim_dict0 = {'UD-COSMOS':  [ 1.30, 3.00],
                    'UD-SXDS':    [-5.80,-4.10],
                    'D-COSMOS':   [ 2.00, 3.75],
                    'D-DEEP2_3':  [-1.90, 1.10],
                    'D-ELAIS_N1': [52.90,56.85]}

    for t_field in gal_field0:
        f_idx = [xx for xx in range(len(tab0)) if gal_field[xx] == t_field]

        fig, ax = plt.subplots()

        n_subsample = 0

        for key in sample_keys:
            s_idx = sub_dict0[key]

            # Get the intersection of f_idx and s_idx
            t_idx = list(set(f_idx) & set(s_idx))
            print t_field, key, len(t_idx)

            if len(t_idx) > 0:
                m0 = 'o' if 'NB0921' in key else '+' if 'NB0816' in key else ''
                c0 = 'red' if 'Ha' in key else 'green' if 'OIII' in key \
                     else 'blue'
                t_name = key.replace('NB0','NB')
                t_name = t_name.replace('OIII_','[OIII] ')
                t_name = t_name.replace('OII_','[OII] ')
                t_name = t_name.replace('Ha_',r'H$\alpha$ ')
                t_name += ' ('+str(len(t_idx))+')'

                # Mod on 01/03/2018
                if (noOII and 'OII_' in key):
                    if silent == False: log.warn('## Will not plot '+key)
                else:
                    n_subsample += 1
                    ax.scatter(ra0[t_idx], dec0[t_idx], s=5, marker=m0,
                               color=c0, linewidth=0.5, edgecolor='none',
                               alpha=0.5, label=t_name)
            #endif
        #endfor

        ax.set_xlabel('Right Ascension (deg)')
        ax.set_ylabel('Declination (deg)')

        ax.set_xlim([RA_lim_dict0[t_field][1], RA_lim_dict0[t_field][0]])
        ax.set_ylim(DE_lim_dict0[t_field])
        ax.annotate(t_field, [0.025,0.975], xycoords='axes fraction',
                    fontsize=12, fontweight='bold', ha='left', va='top')
        ax.minorticks_on()

        # Mod on 01/03/2018
        ncol = 3 if n_subsample % 3 == 0 else 2
        ax.legend(loc='lower center', ncol=ncol, frameon=False, fontsize=10,
                  framealpha=0.9)

        # Overlay DEIMOS FoV | + on 01/03/2018
        if DEIMOS:
            a_coord = [np.average(ra0[f_idx]), np.average(dec0[f_idx])]
            ax = plot_deimos_fov(ax, a_coord, pa=0.0)

        # Overlay Hecto FoV | + on 01/03/2018
        if Hecto:
            a_coord = [np.average(ra0[f_idx]), np.average(dec0[f_idx])]
            ax = plot_hecto_fov(ax, a_coord)

        # Mod on 01/03/2018
        out_pdf = dir0 + 'plots/' + t_field+'_radec.pdf'
        if noOII: out_pdf = out_pdf.replace('.pdf','.noOII.pdf')
        if DEIMOS: out_pdf = out_pdf.replace('.pdf', '.DEIMOS.pdf')
        if Hecto: out_pdf = out_pdf.replace('.pdf', '.Hecto.pdf')
        if silent == False: log.info('## Writing : '+out_pdf)
        fig.savefig(out_pdf, bbox_inches='tight')

    if silent == False: log.info('### End main : '+systime())
#enddef
