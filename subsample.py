"""
subsample
=========

Select NB excess emitters for each redshift sample
"""

from chun_codes import systime

from os.path import exists
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

from astropy.table import Table
from astropy import log

from collections import OrderedDict

from . import read_catalog

def main(tab0=None, field='', dr='pdr1', silent=False, verbose=True):

    '''
    Provide explanation for function here.

    Parameters
    ----------
    tab0 : astropy.table.table.Table
      Astropy Table of HSC-SSP NB excess emitter catalog.
      If not provided, table will be read in based on field and dr inputs

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
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    if type(tab0) == type(None) and field == '':
        log.warn("### [field] input not specified!!!")
        log.warn("### Either 'udeep' or 'deep'")
        log.warn('### Exiting!!!')
        return

    if type(tab0) == type(None):
        tab0 = read_catalog.main(field=field, dr=dr, silent=silent,
                                 verbose=verbose)

    n_gal = len(tab0)

    em_line = [('OII' if xx == 1 else ('OIII' if xx == 3 else 'Ha')) for
               xx in tab0['emission_line']]
    subsam = [a+'_'+b for a,b in zip(em_line, tab0['nb'])]

    subsam0 = list(set(subsam))
    if silent == False:
        log.info('### Here are the subsamples : '+', '.join(subsam0))

    sub_dict0 = OrderedDict()
    for ss in range(len(subsam0)):
        idx = np.array([xx for xx in range(n_gal) if subsam[xx] == subsam0[ss]])
        sub_dict0[subsam0[ss]] = idx

    return sub_dict0
    if silent == False: log.info('### End main : '+systime())
#enddef

