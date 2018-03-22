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

import glob # + on 20/03/2018

from astropy.table import Table
from astropy import log

import matplotlib.patches as mpatches # + on 01/03/2018
from matplotlib.path import Path as Path0 # + on 03/03/2018

# + on 19/03/2018
import shapely.geometry as sg
import descartes

from .. import read_catalog
from .. import subsample
from .. import paths
from . import ds9_mask_overlay # + on 20/03/2018

from .deimos import in_deimos_field, plot_deimos_fov # + on 21/03/2018
from .hecto import in_hecto_field, plot_hecto_fov # + on 21/03/2018
from .bino import plot_bino_fov # + on 22/03/2018
from .subsample_in_pointing import main as ss_in_ptg # + on 21/03/2018

#bbox_props = dict(boxstyle="square,pad=0.15", fc="white", alpha=0.9, ec="none")

def overlay_primus(PRIMUS_tab0, ax):
    '''
    Overlay PRIMUS pointings

    Parameters
    ----------

    Returns
    -------
    ax with Circle patch added

    Notes
    -----
    Created by Chun Ly, 19 March 2018
     - Use shapely and descartes to plot PRIMUS pointing as intersection
       of square (27.2') and circle (radius=14.94')
    Modified by Chun Ly, 21 March 2018
     - Plotting aesthetics (transparent field boundaries)
    '''

    c_ra  = PRIMUS_tab0['RACEN'].data
    c_dec = PRIMUS_tab0['DECCEN'].data

    bsz = 0.453/2.0 # in degree
    for ii in range(len(PRIMUS_tab0)):
        coord = [c_ra[ii], c_dec[ii]]
        circ = sg.Point(coord[0],coord[1]).buffer(0.249) #30.0/60.0)
        box  = sg.box(c_ra[ii]-bsz,c_dec[ii]-bsz,c_ra[ii]+bsz,c_dec[ii]+bsz)
        int0 = circ.intersection(box)
        ax.add_patch(descartes.PolygonPatch(int0, fc='none', ec='purple',
                                            alpha=0.3))
    #endfor

    return ax
#enddef

def main(field='', dr='pdr1', noOII=False, DEIMOS=False, Hecto=False,
         Bino=False, silent=False, verbose=True):

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

    Bino : boolean
      Indicate whether to plot MMT/Bino FoV on PDF plots

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
    Modified by Chun Ly, 4 March 2018
     - Call in_deimos_field()
     - Create table of targets in DEIMOS field, deimos_tab0
    Modified by Chun Ly, 5 March 2018
     - Write table of targets in DEIMOS field to LaTeX files
     - Include RA, Dec, PA in DEIMOS target field table
     - Simplify fld_arr0 cmd for exec
     - Define and pass maskno into plot_deimos_fov()
    Modified by Chun Ly, 7 March 2018
     - Get average RA and Dec for fields
    Modified by Chun Ly, 16 March 2018
     - Overlay Hecto pointings from input coordinate list
     - Hecto sub-code tested. Bug fix: f_idx -> h_idx
    Modified by Chun Ly, 19 March 2018
     - Read in PRIMUS catalog; Call overlay_primus()
    Modified by Chun Ly, 20 March 2018
     - Import and call ds9_mask_overlay() to overlay ds9 regions
    Modified by Chun Ly, 21 March 2018
     - Add Bino boolean keyword and plot Binospec FoV when set
     - Change out_pdf suffix for Bino and Hecto case
     - List comprehensions for mask/config centers
     - Call subsample_in_pointing for DEIMOS fields
     - Call subsample_in_pointing for Hecto fields
    '''

    if silent == False: log.info('### Begin main : '+systime())

    if field == '':
        log.warn("### [field] input not specified!!!")
        log.warn("### Either 'udeep' or 'deep'")
        log.warn('### Exiting!!!')
        return

    dir0 = paths.gdrive() # Mod on 02/03/2018

    # + on 03/03/2018
    if DEIMOS or Hecto or Bino:
        field_coord_file = dir0 + 'field_coordinates.txt'
        if silent == False: log.info('### Writing : '+field_coord_file)
        ptg_tab = asc.read(field_coord_file)

    # + on 19/03/2018
    PRIMUS_cat_file = dir0 + 'catalogs/primus/PRIMUS_mask_centers.txt'
    PRIMUS_tab0     = asc.read(PRIMUS_cat_file)

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

        # + on 07/03/2018
        RA_avg = np.average(ra0[f_idx])
        DE_avg = np.average(dec0[f_idx])
        print 'RA/DE avg : ', t_field, RA_avg, DE_avg

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

        # Overlay PRIMUS pointings | + on 19/03/2018
        ax = overlay_primus(PRIMUS_tab0, ax)

        # Overlay ds9 masked regions | + on 20/03/2018
        tmp_field = t_field.split('-')[-1].replace('_','-')
        prefix = 'mask_dr1_s15b_'+field+'_NB0921_'+tmp_field+'*.reg'
        ds9_file = glob.glob(dir0+'catalogs/masks_Hayashi+17/'+prefix)[0]
        if silent == False: log.info('## ds9_file : '+ds9_file)
        ds9_mask_overlay.main(ax, ds9_file, color='black', alpha=0.15)

        # Overlay DEIMOS FoV | + on 01/03/2018, Mod on 03/03/2018
        if DEIMOS:
            d_idx = [xx for xx in range(len(ptg_tab)) if
                     ((ptg_tab['Instr'][xx] == 'DEIMOS') and
                      (ptg_tab['Field'][xx] == t_field))]
            if len(d_idx) > 0:
                t_tab = ptg_tab[d_idx]
                pa = t_tab['PA'].data
                a_coord = [[ra,de] for ra,de in zip(t_tab['RA'], t_tab['Dec'])]

                # Mod on 05/03/2018
                maskno = [mname.replace(t_field+'-D','') for
                          mname in t_tab['MaskName']]
                ax, deimos_verts0 = plot_deimos_fov(ax, a_coord, maskno, pa=pa)

                # + on 04/03/2018
                deimos_fld_idx = in_deimos_field(tab0, deimos_verts0,
                                                 silent=silent, verbose=verbose)

                # Get subsample sizes in each DEIMOS pointing | + on 04/03/2018
                # Mod on 21/03/2018
                deimos_outfile = dir0+'catalogs/'+t_field+'_deimos.tex'
                ss_in_ptg(sub_dict0, deimos_fld_idx, t_tab, 'DEIMOS',
                          deimos_outfile, silent=silent, verbose=verbose)
            #endif
        #endif

        # Overlay Hecto FoV | + on 01/03/2018, Mod on 16/03/2018
        if Hecto:
            h_idx = [xx for xx in range(len(ptg_tab)) if
                     ((ptg_tab['Instr'][xx] == 'Hecto') and
                      (ptg_tab['Field'][xx] == t_field))]
            if len(h_idx) > 0:
                t_tab = ptg_tab[h_idx]
                a_coord = [[ra,de] for ra,de in zip(t_tab['RA'], t_tab['Dec'])]

                configno = [cname.replace(t_field+'-H','') for
                            cname in t_tab['MaskName'].data]
                ax = plot_hecto_fov(ax, a_coord, configno)

                hecto_fld_idx = in_hecto_field(tab0, a_coord, silent=silent,
                                               verbose=verbose)

                # Get subsample sizes in each Hecto pointing | Mod on 21/03/2018
                hecto_outfile = dir0+'catalogs/'+t_field+'_hecto.tex'
                ss_in_ptg(sub_dict0, hecto_fld_idx, t_tab, 'Hecto',
                          hecto_outfile, silent=silent, verbose=verbose)
            #endif
        #endif

        # Overlay Binospec FoV | + on 21/03/2018
        if Bino:
            b_idx = [xx for xx in range(len(ptg_tab)) if
                     ((ptg_tab['Instr'][xx] == 'Bino') and
                      (ptg_tab['Field'][xx] == t_field))]
            if len(b_idx) > 0:
                t_tab = ptg_tab[b_idx]
                pa = t_tab['PA'].data
                a_coord = [[ra,de] for ra,de in zip(t_tab['RA'], t_tab['Dec'])]

                maskno = [mname.replace(t_field+'-B','') for
                          mname in t_tab['MaskName'].data]
                ax, deimos_verts0 = plot_bino_fov(ax, a_coord, maskno, pa=pa)
            #endif
        #endif

        # Mod on 01/03/2018
        out_pdf = dir0 + 'plots/' + t_field+'_radec.pdf'
        if noOII: out_pdf = out_pdf.replace('.pdf','.noOII.pdf')
        if DEIMOS: out_pdf = out_pdf.replace('.pdf', '.DEIMOS.pdf')
        # Mod on 21/03/2018
        if Hecto and Bino:
            out_pdf = out_pdf.replace('.pdf', '.MMT.pdf')
        else:
            if Hecto: out_pdf = out_pdf.replace('.pdf', '.Hecto.pdf')
            if Bino: out_pdf = out_pdf.replace('.pdf', '.Bino.pdf')

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
