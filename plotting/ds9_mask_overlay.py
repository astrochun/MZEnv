"""
ds9_mask_overlay
================

Overlay ds9 regions as matplotlib patches
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

def main(ax0, maskregfile, color='black', alpha=1.0, zorder=0):
    '''
    Main function of ds9_mask_overlay

    Parameters
    ----------
    ax0 : matplotlib.axes._subplots.AxesSubplot
      Matplotlib axes from plt.subplots()

    maskregfile : str
      Full path filename for ds9 region file to overlay

    color : str
      Color to use for shading. Default: Black

    alpha : float
      Transparency value. (0 = transparent, 1 = opaque). Default: 1

    zorder : int
      Matplotlib zorder for overlaying levels

    Returns
    -------
    ax0 with ds9 patches added

    Notes
    -----
    Created by Chun Ly, 20 March 2018
     - This was in fact provided by Masao Hayashi with minor modifications
    '''
    
    f = open(maskregfile,'r')
    mask = []
    for line in f:
        if 'box' == line[0:3] or 'circle' == line[0:6]:
            tmpmask = []
            tmp1 = line.split('(')
            tmpmask.append(tmp1[0])
            tmp2 = tmp1[1].split(')')[0].split(',')
            tmpmask.extend([tmp2[0],tmp2[1]])
            del tmp2[0]
            del tmp2[0]
            tmpmask.extend([','.join(tmp2)])
            mask.append(tuple(tmpmask))
        #endif
    #endfor
    mask = np.array(mask,dtype=[('shape','a8'),('ra2000','f8'),
                                ('decl2000','f8'), ('parameters','a128')])
    f.close()

    for tmp in mask:
        ra    = float(tmp['ra2000'])
        dec   = float(tmp['decl2000'])
        if tmp['shape'] == 'box':
            tmp_param = tmp['parameters'].split(',')
            dx    = float(tmp_param[0][:-1])/3600.0 / np.cos(np.radians(dec))
            dy    = float(tmp_param[1][:-1])/3600.0
            theta = float(tmp_param[2])
            rad0  = np.radians(-theta)
            # ra2,dec2 is at lower left corner
            #        dec
            #         ^
            #         |
            #   ra <--    thus, -theta
            #
            ra2  = ra - ( 0.5 * dx * np.cos(rad0)  - 0.5 * dy * np.sin(rad0) ) 
            dec2 = dec - ( 0.5 * dx * np.sin(rad0) + 0.5 * dy * np.cos(rad0) )
            rectangle = plt.Rectangle([ra2,dec2], dx, dy, -theta, ec='none',
                                      fc=color, alpha=alpha)
            ax0.add_patch(rectangle)
        elif tmp['shape'] == 'circle':
            r = float(tmp['parameters'][:-1])/3600.0
            circle = Ellipse(xy=(ra,dec), width=r/np.cos(np.radians(dec))*2.0,
                             height=r*2.0, angle=0.0, ec='none', fc=color,
                             alpha=alpha)
            ax0.add_patch(circle)
    #endfor

    return ax0
#enddef
