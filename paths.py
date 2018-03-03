"""
paths
=====

Defines main paths for each user
"""

from getpass import getuser
from astropy import log

def gdrive(silent=True, verbose=False):

    '''
    Google Drive path for MZEnv project

    Parameters
    ----------

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 2 March 2018
    '''
    
    if silent == False: log.info('### Begin gdrive')

    if getuser() == 'cly':
        dir0 = '/Users/cly/Google Drive/MZEnv/'

    if silent == False: log.info('### End gdrive')
    return dir0
#enddef

