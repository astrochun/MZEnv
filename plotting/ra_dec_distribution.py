"""
ra_dec_distribution
===================

Plot spatial distribution of HSC-SSP NB excess emitter samples
"""

import sys, os

from chun_codes import systime

from os.path import exists
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import log

from .. import read_catalog

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
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    if field == '':
        log.warn("### [field] input not specified!!!")
        log.warn("### Either 'udeep' or 'deep'")
        log.warn('### Exiting!!!')
        return

    tab0 = read_catalog.main(field=field, dr=dr, silent=silent, verbose=verbose)

    if silent == False: log.info('### End main : '+systime())
#enddef

