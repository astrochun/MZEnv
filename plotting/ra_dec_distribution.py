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
from matplotlib.path import Path as Path0 # + on 03/03/2018

from .. import read_catalog
from .. import subsample
from .. import paths

#bbox_props = dict(boxstyle="square,pad=0.15", fc="white", alpha=0.9, ec="none")

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


def plot_deimos_fov(ax, coord, pa=0.0):
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
    '''

    dx0 =  5.0 / 60.0 # in deg
    dy0 = 16.7 / 60.0 # in deg

    xbox0 = np.array([ 1, 1,-1,-1])*dx0/2.0
    ybox0 = np.array([-1, 1, 1,-1])*dy0/2.0

    n_ptgs = len(coord) # + on 03/03/2018

    verts0 = [] # + on 03/03/2018

    # Mod on 03/03/2018
    for cc in range(n_ptgs):
        rad = -1*pa[cc]*np.pi/180.0 #Sign flip

        xprime =  xbox0*np.cos(rad) - ybox0*np.sin(rad) #Sign fix
        yprime =  xbox0*np.sin(rad) + ybox0*np.cos(rad) #Sign fix

        xbox = coord[cc][0] + xprime
        ybox = coord[cc][1] + yprime

        xy_low = [xbox[3], ybox[3]]
        ax.scatter([coord[cc][0]], [coord[cc][1]])
        rect = mpatches.Rectangle(xy_low, dx0, dy0, angle=-1*pa[cc], alpha=0.33,
                                  ec="r", color="r") # Sign flip

        # Must be done before ax.add_patch call | + on 03/03/2018
        verts0.append(rect.get_verts())
        ax.add_patch(rect)
    #endfor

    return ax, verts0 # Mod on 03/03/2018
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
     - Aesthetic changes (unfilled circle)
    '''

    radius = 0.5 # in deg

    circ = mpatches.Circle(coord, radius=radius, alpha=0.5, ec="b",
                           color="none", linewidth=1.5)
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
     - Update call to subsample.main(); simplify code
    Modified by Chun Ly, 2 March 2018
     - Call paths module to get path
    Modified by Chun Ly, 3 March 2018
     - Read in field coordinate file for specific FoV overlay
     - Generate DEIMOS coordinate list, pass to plot_deimos_fov()
     - Bug fix: Call plot_deimos_fov() outside of for loop
     - List comprehension: require DEIMOS pointing inside field
     - If statement before calling plot_deimos_fov()
     - Get vertices from plot_deimos_fov()
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    if field == '':
        log.warn("### [field] input not specified!!!")
        log.warn("### Either 'udeep' or 'deep'")
        log.warn('### Exiting!!!')
        return

    dir0 = paths.gdrive() # Mod on 02/03/2018

    # + on 03/03/2018
    if DEIMOS or Hecto:
        field_coord_file = dir0 + 'field_coordinates.txt'
        if silent == False: log.info('### Writing : '+field_coord_file)
        ptg_tab = asc.read(field_coord_file)

    tab0 = read_catalog.main(field=field, dr=dr, silent=silent, verbose=verbose)
    ra0, dec0 = tab0['ra'], tab0['dec']

    # Mod on 01/03/2018
    sub_dict0, gal_dict0 = subsample.main(tab0=tab0, field=field, dr=dr)

    gal_field0 = gal_dict0.keys() # Unique galaxy field list

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
        f_idx = gal_dict0[t_field] # Mod on 01/03/2018

        fig, ax = plt.subplots()

        n_subsample = 0

        for key in sub_dict0.keys(): # Mod on 01/03/2018
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

        # Overlay DEIMOS FoV | + on 01/03/2018, Mod on 03/03/2018
        if DEIMOS:
            d_idx = [xx for xx in range(len(ptg_tab)) if
                     ((ptg_tab['Instr'][xx] == 'DEIMOS') and
                      (ptg_tab['Field'][xx] == t_field))]
            if len(d_idx) > 0:
                pa, a_coord = [], []
                for cc in range(len(d_idx)):
                    t_tab = ptg_tab[d_idx[cc]]
                    a_coord.append([t_tab['RA'], t_tab['Dec']])
                    pa.append(t_tab['PA'])

                ax, deimos_verts0 = plot_deimos_fov(ax, a_coord, pa=pa)

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


def run_18b():

    '''
    Call main() for specific plots to generate sets

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Created by Chun Ly, 1 March 2018
    '''

    for depth in ['udeep','deep']:
        main(depth, dr='pdr1', DEIMOS=True)
        main(depth, dr='pdr1', Hecto=True)
#enddef
