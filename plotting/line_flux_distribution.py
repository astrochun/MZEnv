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

    out_pdf = dir0+'plots/'+field+'_line_flux.pdf'
    if silent == False: log.info('Output PDF : '+out_pdf)
    pp = PdfPages(out_pdf)

    gal_field0 = gal_dict0.keys() # Unique galaxy field list

    emline = [r'H$\alpha$',r'[OIII]$\lambda$5007',r'[OII]$\lambda$3727']
    ctype = ['r', 'g', 'b'] # for Ha, [OIII], [OII]

    for t_field in gal_field0:
        f_idx = gal_dict0[t_field]

        fig, ax_arr = plt.subplots(nrows=3, ncols=2)
        fig.suptitle(t_field, fontsize=14)

        for key in sub_dict0.keys():
            s_idx = sub_dict0[key]

            # Get the intersection of f_idx and s_idx
            t_idx = list(set(f_idx) & set(s_idx))
            print t_field, key, len(t_idx)

            logFlux = tab0['line_flux_observed']
            flux_min = np.min(logFlux[s_idx])
            flux_max = np.max(logFlux[s_idx])

            binsize = 0.1
            f_start = flux_min - (flux_min - np.floor(flux_min)) % binsize
            f_end   = flux_max + (flux_max - np.floor(flux_max)) % binsize
            bins    = np.arange(f_start,f_end, binsize)
            if len(t_idx) > 0:
                Flux = logFlux[t_idx]

                if 'Ha_'   in key: row = 0
                if 'OIII_' in key: row = 1
                if 'OII_'  in key: row = 2

                if 'NB0816' in key: col = 0
                if 'NB0921' in key: col = 1

                t_ax = ax_arr[row,col]
                if row == 0:
                    col_title = key.replace('Ha_','').replace('NB0','NB')
                    t_ax.set_title(col_title, fontsize=12)

                t_ax.hist(Flux, bins=bins, color=ctype[row], edgecolor='none',
                          alpha=0.5, histtype='stepfilled', align='mid')
                if row < 2:
                    t_ax.xaxis.set_ticklabels([])
                else:
                    t_ax.set_xlabel(r'$\log({\rm Flux}_{\rm NB})$')

                t_ax.annotate(emline[row], [0.975,0.95], va='top', ha='right',
                              xycoords='axes fraction', fontsize=10)
                t_ax.minorticks_on()
            #endif
        #endfor

        subplots_adjust(left=0.025, bottom=0.025, top=0.90, right=0.975,
                        wspace=0.13, hspace=0.02)
        fig.set_size_inches(8.0, 10)
        fig.savefig(pp, format='pdf', bbox_inches='tight')
        plt.close()
        fig.clear()

    #endfor

    if silent == False: log.info('Writing : '+out_pdf)
    pp.close()
    
    if silent == False: log.info('### End main : '+systime())
#enddef

