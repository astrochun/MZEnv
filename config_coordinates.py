from astropy.io import ascii as asc

from astropy.coordinates import SkyCoord
import astropy.units as u

from astropy.table import Column

from . import paths

def main():
    dir0 = paths.gdrive() # Mod on 02/03/2018

    tab0 = asc.read(dir0+'field_coordinates.txt')

    c0 = SkyCoord(ra=tab0['RA'], dec=tab0['Dec'], unit=(u.deg, u.deg))
    hmsdms = c0.to_string('hmsdms')
    mmt_idx = [xx for xx in range(len(tab0)) if tab0['Instr'][xx] == 'Bino'
               or tab0['Instr'][xx] == 'Hecto']

    hms = [str0.split(' ')[0].replace('h',':').replace('m',':').replace('s','')
           for str0 in hmsdms]
    dms = [str0.split(' ')[1].replace('d',':').replace('m',':').replace('s','')
           for str0 in hmsdms]

    ch = Column(hms, name='HMS')
    cd = Column(dms, name='DMS')

    tab0.add_columns([ch, cd])
    tab0[mmt_idx].pprint(max_lines=-1)

