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

from .. import read_catalog
from .. import subsample

#bbox_props = dict(boxstyle="square,pad=0.15", fc="white", alpha=0.9, ec="none")

def main(field='', dr='pdr1', silent=False, verbose=True):

    '''
    Main function to plot RA and Dec for each sub-sample of NB excess emitters

    Parameters
    ----------
    field : str
      Name of field to read in. Either 'udeep' or 'deep'

    dr : str
      Data release name. Default: 'pdr1'

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

    for t_field in gal_field0:
        f_idx = [xx for xx in range(len(tab0)) if gal_field[xx] == t_field]

        fig, ax = plt.subplots()

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
                ax.scatter(ra0[t_idx], dec0[t_idx], s=5, marker=m0, color=c0,
                           linewidth=0.5, edgecolor='none', alpha=0.5,
                           label=t_name)
            #endif
        #endfor

        ax.set_xlabel('Right Ascension (deg)')
        ax.set_ylabel('Declination (deg)')

        ax.annotate(t_field, [0.05,0.95], xycoords='axes fraction',
                    fontsize=12, ha='left', va='top')
        ax.minorticks_on()
        ax.legend(loc='lower right', frameon=False, fontsize=10, framealpha=0.9)

        out_pdf = dir0 + 'plots/' + t_field+'_radec.pdf'
        fig.savefig(out_pdf, bbox_inches='tight')

    if silent == False: log.info('### End main : '+systime())
#enddef

