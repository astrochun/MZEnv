"""
line_flux_distribution
======================

Generate emission-line fluxes based on NB excess emission
"""

import sys, os

from chun_codes import systime

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pylab import subplots_adjust

from astropy import log

from .. import read_catalog
from .. import subsample
from .. import paths

deimos_limit = np.log10(5e-18) # 5-sigma for 1.5 hr

def sens_overlay(t_ax, limit, text0, ymin=0.90, ymax=0.975, color='k'):
    '''
    Draw vertical lines for sensitivity and label them

    Parameters
    ----------
    t_ax : matplotlib.axes._subplots.AxesSubplot
      Matplotlib axes from plt.subplots()

    limit : flux limit value in log units

    text0 : str
      String for annotation

    ymin : float
      Minimum value to draw vertical line. Default: 0.90

    ymax : float
      Maximum value to draw vertical line. Default: 0.975

    color : str
      Color of line. Default: black

    Returns
    -------
      t_ax with vertical line and text annotation

    Notes
    -----
    Created by Chun Ly, 13 March 2018
    '''

    t_ax.axvline(limit, ymin=ymin, ymax=ymax, color=color)
    ax_ymin, ax_ymax = t_ax.get_ylim()
    t_ax.text(limit+0.05, ax_ymax*ymax, text0, fontsize=8, ha='left', va='top')

    return t_ax
#enddef

def main(field='', dr='pdr1', DEIMOS=False, Hecto=False, silent=False,
         verbose=True):

    '''
    Main function to generate histogram plots for various emission lines

    Parameters
    ----------
    field : str
      Name of field to read in. Either 'udeep' or 'deep'

    dr : str
      Data release name. Default: 'pdr1'

    DEIMOS : boolean
      Indicate whether to plot Keck/DEIMOS sensitivity

    Hecto : boolean
      Indicate whether to plot MMT/Hecto sensitivity

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 13 March 2018
     - Plotting aesthetics: Same binning and x limits for all em lines
     - Plotting aesthetics: Remove subplots if no data
     - Plotting aesthetics: White space improvements
     - Plot DEIMOS sensitivity when specified
     - Call sens_overlay() function for Hb line
     - Handle Hb line overlay for Ha and [OIII] emitters
     - Add [OIII]5007 sensitivity for Ha emitters
     - Add [OIII]4363 sensitivity for Ha and [OIII] emitters
     - Add [OII] sensitivity for Ha, [OIII], and [OII] emitters
    Modified by Chun Ly, 14 March 2018
     - Determine fraction of sources detected for various emission lines
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    if field == '':
        log.warn("### [field] input not specified!!!")
        log.warn("### Either 'udeep' or 'deep'")
        log.warn('### Exiting!!!')
        return

    dir0 = paths.gdrive() # Mod on 02/03/2018

    tab0 = read_catalog.main(field=field, dr=dr, silent=silent, verbose=verbose)

    sub_dict0, gal_dict0 = subsample.main(tab0=tab0, field=field, dr=dr)

    gal_field0 = gal_dict0.keys() # Unique galaxy field list

    # + on 14/03/2018
    det_Hb    = np.zeros( (len(gal_field0), len(sub_dict0.keys())) )
    det_Hg    = np.zeros( (len(gal_field0), len(sub_dict0.keys())) )
    det_OIIIa = np.zeros( (len(gal_field0), len(sub_dict0.keys())) )

    out_pdf = dir0+'plots/'+field+'_line_flux.pdf'
    if silent == False: log.info('Output PDF : '+out_pdf)
    pp = PdfPages(out_pdf)

    gal_field0 = gal_dict0.keys() # Unique galaxy field list

    emline = [r'H$\alpha$',r'[OIII]$\lambda$5007',r'[OII]$\lambda$3727']
    ctype = ['r', 'g', 'b'] # for Ha, [OIII], [OII]

    ff = 0
    for t_field in gal_field0:
        f_idx = gal_dict0[t_field]

        logFlux = tab0['line_flux_observed']
        flux_min = np.min(logFlux[f_idx])
        flux_max = np.max(logFlux[f_idx])

        binsize = 0.1
        f_start = flux_min - (flux_min - np.floor(flux_min)) % binsize - 0.10
        f_end   = flux_max + (flux_max - np.floor(flux_max)) % binsize + 0.05
        bins0   = np.arange(f_start,f_end, binsize)

        fig, ax_arr = plt.subplots(nrows=3, ncols=2)
        fig.suptitle(t_field, fontsize=14)

        ss = 0 # Count for indexing | + on 14/03/2018
        for key in sub_dict0.keys():
            s_idx = sub_dict0[key]

            # Get the intersection of f_idx and s_idx
            t_idx = list(set(f_idx) & set(s_idx))
            print t_field, key, len(t_idx)

            if 'Ha_'   in key: row = 0
            if 'OIII_' in key: row = 1
            if 'OII_'  in key: row = 2

            if 'NB0816' in key: col = 0
            if 'NB0921' in key: col = 1

            t_ax = ax_arr[row,col]

            if len(t_idx) > 0:
                Flux = logFlux[t_idx]

                if row == 0:
                    col_title = key.replace('Ha_','').replace('NB0','NB')
                    t_ax.set_title(col_title, fontsize=12)

                t_ax.hist(Flux, bins=bins0, color=ctype[row], edgecolor='none',
                          alpha=0.5, histtype='stepfilled', align='mid')
                t_ax.set_xlim([f_start, f_end])
                if row < 2:
                    t_ax.xaxis.set_ticklabels([])
                else:
                    t_ax.set_xlabel(r'$\log({\rm Flux}_{\rm NB})$')

                t_ax.annotate(emline[row], [0.975,0.95], va='top', ha='right',
                              xycoords='axes fraction', fontsize=10)
                t_ax.minorticks_on()

                if DEIMOS:
                    if row == 0:
                        # This is for A(Ha) = 1 and 10-sigma limit
                        Hb_lim = deimos_limit + np.log10(4.22 * 10/5.)
                        text0 = r'H$\beta$ S/N=10'+'\n'+r'A(H$\alpha$)=1'
                    if row == 1:
                        # This is for [OIII]/Hb = 5 and 10-sigma limit
                        Hb_lim = deimos_limit + np.log10(5 * 10/5.)
                        text0 = r'H$\beta$ S/N=10'+'\n'+r'[OIII]/H$\beta$=5'

                    if row <= 1:
                        sens_overlay(t_ax, Hb_lim, text0, ymin=0.90, ymax=0.975,
                                     color='r')

                        # + on 14/03/2018
                        Hb_lim5 = Hb_lim - np.log10(2.0) # This is for 5-sigma
                        Hb_idx = [xx for xx in range(len(Flux)) if
                                  Flux[xx] >= Hb_lim5]
                        det_Hb[ff,ss] = np.float(len(Hb_idx))/len(Flux) * 100.0

                    if row == 0:
                        # OIII/Ha = 1 and 100-sigma limit
                        OIII_lim = deimos_limit + np.log10(100/5.)
                        text0 = r'[OIII] S/N=100'+'\n'+r'[OIII]/H$\alpha$=1'
                        sens_overlay(t_ax, OIII_lim, text0, ymin=0.80,
                                     ymax=0.875, color='g')

                    if row <=1:
                        # OIII/Ha = 1 and 4363/5007=0.01, 3-sigma limit
                        OIIIa_lim = deimos_limit + np.log10(3./5/0.01)
                        text0 = r'[OIII]$\lambda$4363 S/N=3'+'\n'+\
                                r'[OIII]$\lambda$4363/[OIII]$\lambda$5007=0.01'
                        sens_overlay(t_ax, OIIIa_lim, text0, ymin=0.70,
                                     ymax=0.775, color='g')

                        # + on 14/03/2018
                        OIIIa_idx = [xx for xx in range(len(Flux)) if
                                     Flux[xx] >= OIIIa_lim]
                        det_OIIIa[ff,ss] = np.float(len(OIIIa_idx))/len(Flux)*100.0

                    # [OII] fluxes
                    if row == 0:
                        OII_lim = deimos_limit + np.log10(2.0 * 10/5.)
                        text0 = r'[OII] S/N=10'+'\n'+r'[OII]/H$\alpha$=0.5'
                    if row == 1:
                        OII_lim = deimos_limit + np.log10(5.0 * 10/5.)
                        text0 = r'[OII] S/N=10'+'\n'+r'[OII]/[OIII]=0.2'
                    if row == 2:
                        OII_lim = deimos_limit + np.log10(50/5.)
                        text0 = r'[OII] S/N=50'
                    sens_overlay(t_ax, OII_lim, text0, ymin=0.80,
                                 ymax=0.875, color='b')
                #endif
            else:
                t_ax.axis('off')

            ss += 1 # + on 14/03/2018
        #endfor

        subplots_adjust(left=0.025, bottom=0.025, top=0.95, right=0.975,
                        wspace=0.13, hspace=0.04)
        fig.set_size_inches(8.0, 10)
        fig.savefig(pp, format='pdf', bbox_inches='tight')
        plt.close()
        fig.clear()
        ff += 1 # + on 14/03/2018
    #endfor

    if silent == False: log.info('Writing : '+out_pdf)
    pp.close()
    
    if silent == False: log.info('### End main : '+systime())
#enddef

