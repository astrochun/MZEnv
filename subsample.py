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
from getpass import getuser

from . import read_catalog
from . import paths

def galaxy_field(tab0, field, silent=False, verbose=True):
    '''
    Determine galaxy field (e.g., COSMOS, SXDS, etc.) based on RA

    Parameters
    ----------
    tab0 : astropy.table.table.Table
      Astropy Table of HSC-SSP NB excess emitter catalog.

    field : str
      Name of field to read in. Either 'udeep' or 'deep'

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
    g_field0 : numpy.array
      String array containing galaxy field name (e.g., 'UD-COSMOS')

    Notes
    -----
    Created by Chun Ly, 28 February 2018

    Modified by Chun Ly, 1 March 2018
     - Define galaxy dictionary and return instead
     - Bug fix: typo
    '''

    ra  = tab0['ra'] / 15.0 # in hours
    dec = tab0['dec']

    t_field    = np.chararray(len(tab0), itemsize=10)
    t_field[:] = 'N/A'

    sxds = [xx for xx in range(len(tab0)) if
            (ra[xx] >= 2.0 and ra[xx] <= 3.0)]
    if len(sxds) > 0:
        t_field[sxds] = 'SXDS'

    cosmos = [xx for xx in range(len(tab0)) if
              (ra[xx] >= 9.75 and ra[xx] <= 10.25)]
    if len(cosmos) > 0:
        t_field[cosmos] = 'COSMOS'

    elais = [xx for xx in range(len(tab0)) if
             (ra[xx] >= 15.75 and ra[xx] <= 16.5)]
    if len(elais) > 0:
        t_field[elais] = 'ELAIS_N1'

    deep2 = [xx for xx in range(len(tab0)) if
             (ra[xx] >= 23.25 and ra[xx] <= 23.75)]
    if len(deep2) > 0:
        t_field[deep2] = 'DEEP2_3'


    if field == 'udeep': prefix = 'UD'
    if field == 'deep': prefix = 'D'

    g_field = np.array([prefix+'-'+a for a in t_field])

    g_field0 = list(set(g_field))

    g_dict0 = OrderedDict()
    for gal in g_field0:
        idx = np.array([xx for xx in range(len(tab0)) if g_field[xx] == gal])
        g_dict0[gal] = idx

    return g_dict0 #, g_field0
#enddef

def summary_table(sub_dict0, gal_dict0, silent=False, verbose=True):
    '''
    Generate Astropy table of number of subsamples for each galaxy field

    Parameters
    ----------
    sub_dict0 : dict
      Dictionary containing indices for each NB excess emitter subsample.
      From main()

    gal_dict0 : dict
      Dictionary containing indices for each galaxy field name
      [from galaxy_field()]

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
    tab0 : astropy.table.table.Table
      Table of number of excess emitters for each subsample and galaxy field

    Notes
    -----
    Created by Chun Ly, 1 March 2018
     - Change from gal_field to gal_dict0 input
    '''

    gal_field0 = gal_dict0.keys()
    n_fields   = len(gal_field0)

    fieldname   = np.array(gal_field0)
    n_HaNB816   = np.zeros(n_fields, dtype=np.int16)
    n_OIIINB816 = np.zeros(n_fields, dtype=np.int16)
    n_OIINB816  = np.zeros(n_fields, dtype=np.int16)

    n_HaNB921   = np.zeros(n_fields, dtype=np.int16)
    n_OIIINB921 = np.zeros(n_fields, dtype=np.int16)
    n_OIINB921  = np.zeros(n_fields, dtype=np.int16)

    for ff in range(n_fields):
        f_idx = gal_dict0[gal_field0[ff]]

        for key in sub_dict0.keys():
            s_idx = sub_dict0[key]

            # Get the intersection of f_idx and s_idx
            t_idx = list(set(f_idx) & set(s_idx))

            cmd = 'n_%s[%i] = %i ' % (key.replace('_NB0','NB'), ff, len(t_idx))
            exec(cmd)
        #endfor
    #endfor

    arr0 = [fieldname, n_HaNB816, n_OIIINB816, n_OIINB816,
            n_HaNB921, n_OIIINB921, n_OIINB921]
    names0 = ['Field', 'HaNB816', 'OIIINB816', 'OIINB816',
              'HaNB921', 'OIIINB921', 'OIINB921']
    tab0 = Table(arr0, names=names0)

    if silent == False: tab0.pprint()
    return tab0
#enddef

def main(tab0=None, field='', dr='pdr1', silent=False, verbose=True):

    '''
    Main function for HSC NB excess emitter subsample definition

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
    sub_dict0 : dict
      Dictionary containing indices for each NB excess emitter subsample.
      From main()

    gal_field : list
      List of string containing the field name [from galaxy_field()]

    Notes
    -----
    Created by Chun Ly, 28 February 2018
     - Call galaxy_field to get field names for each target
     - Always require [field] (for galaxy_field() input)
    Modified by Chun Ly, 1 March 2018
     - Update documentation
     - Update call to galaxy_field() and return
    Modified by Chun Ly, 2 March 2018
     - Call summary_table()
     - Call paths module to get path
     - Write ASCII-formatted subsample summary table
    Modified by Chun Ly, 6 March 2018
     - Sort by name so Ha is listed before OIII and then OII
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    if field == '':
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
    subsam0.sort() # + on 06/03/2018

    if silent == False:
        log.info('### Here are the subsamples : '+', '.join(subsam0))

    sub_dict0 = OrderedDict()
    for ss in range(len(subsam0)):
        idx = np.array([xx for xx in range(n_gal) if subsam[xx] == subsam0[ss]])
        sub_dict0[subsam0[ss]] = idx

    gal_dict0 = galaxy_field(tab0, field)

    dir0 = paths.gdrive() # Mod on 02/03/2018

    s_tab0 = summary_table(sub_dict0, gal_dict0, silent=silent, verbose=verbose)
    out_summary_tab_file = dir0 + 'catalogs/'+field+'_sample.tbl'
    if silent == False: log.info('### Writing : '+out_summary_tab_file)
    s_tab0.write(out_summary_tab_file, format='ascii.fixed_width_two_line')

    if silent == False: log.info('### End main : '+systime())

    return sub_dict0, gal_dict0
#enddef

