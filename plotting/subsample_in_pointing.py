"""
subsample_in_pointing
=====================

Get sub-sample numbers in each targeted pointing
"""

from chun_codes import systime

import numpy as np

from astropy.table import Table
from astropy import log

def main(sub_dict0, fld_idx, fld_tab, Inst, tab_outfile, silent=False,
         verbose=True):
    '''
    Provide explanation for function here.

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
    Created by Chun Ly, 21 March 2018
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    for key in sub_dict0.keys():
        s_idx = sub_dict0[key]
        
        cmd1 = 'n_fld_%s = np.zeros(len(fld_idx), dtype=np.int)' % key
        exec(cmd1)
        for ff,in_field in enumerate(fld_idx):
            in_idx = list(set(in_field.tolist()) & set(s_idx))
            cmd2 = 'n_fld_%s[ff] = len(in_idx)' % key
            exec(cmd2)
        #endfor
    #endfor
    t_cols = ['n_fld_'+ aa for aa in sub_dict0.keys()]

    fld_arr0 = [fld_tab['MaskName'].data, fld_tab['RA'].data,
                fld_tab['Dec'].data, fld_tab['PA'].data]
    names0 = ['MaskName', 'RA', 'Dec', 'PA']

    cmd3 = "fld_arr0 += ["+', '.join(t_cols)+']'
    exec(cmd3)
    names0 += [val.replace('NB0','NB') for val in sub_dict0.keys()]
    inptg_tab0 = Table(fld_arr0, names=names0)

    if Inst == 'Hecto': fld_arr0.remove_column('PA')
    
    if silent == False: inptg_tab0.pprint(max_lines=-1)
    if silent == False: log.info('### Writing : '+tab_outfile)
    inptg_tab0.write(tab_outfile, format='ascii.latex')

    if silent == False: log.info('### End main : '+systime())
#enddef

