"""
read_catalog
============

Read in HSC SSP catalogs of NB excess emitters
"""

from getpass import getuser
from os.path import exists

from chun_codes import systime

import astropy.io.ascii as asc
import astropy.io.fits as fits

import numpy as np

import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import log

from . import paths

def main(field='', dr='pdr1', silent=False, verbose=True):

    '''
    Main function to read in HSC SSP catalogs of NB excess emitters

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
    Created by Chun Ly, 27 February 2018

    Modified by Chun Ly, 2 March 2018
     - Call paths module to get path
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    if field == '':
        log.warn("### [field] input not specified!!!")
        log.warn("### Either 'udeep' or 'deep'")
        log.warn('### Exiting!!!')
        return

    dir0 = paths.gdrive() # Mod on 02/03/2018

    infile = dir0+'catalogs/%s_%s_nb_forced_coord_phot.csv.gz' % (dr, field)
    if silent == False: log.info('### Reading : '+infile)
    tab0 = asc.read(infile, format='csv')

    if silent == False: log.info('### End main : '+systime())

    return tab0
#enddef

