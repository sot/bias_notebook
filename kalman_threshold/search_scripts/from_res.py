import numpy as np
from astropy.table import Table
import mica.quaternion
from Ska.quatutil import radec2eci
import agasc
from kadi import events
from Ska.engarchive import fetch
from Ska.Matplotlib import plot_cxctime
from Ska.Numpy import smooth
import mica.starcheck
from Chandra.Time import DateTime
import matplotlib.pyplot as plt


msids = ['AOACASEQ', 'AOPCADMD', 'AOKALSTR',  'COTLRDSF', 'COBSRQID',
         'AOATTQT1', 'AOATTQT2', 'AOATTQT3', 'AOATTQT4']

res_msids = ['AORESY0','AORESY1_1', 'AORESY2_1', 'AORESY3', 'AORESY4',
             'AORESY5_1','AORESY6_1', 'AORESY7',
             'AORESZ0', 'AORESZ1_1',
             'AORESZ2_1','AORESZ3', 'AORESZ4', 'AORESZ5_1', 'AORESZ6_1',
             'AORESZ7']
slot_cols = ['AOACFID', 'AOACMAG', 'AOACYAN', 'AOACZAN', 'AOACFCT', 'AOIMAGE',
             'AOACIDP', 'AOACIIR', 'AOACIMS', 'AOACISP']
slots = range(0, 8)
slot_msids = ["{}{}".format(field, slot) for field in slot_cols
                  for slot in slots]


def consecutive(data, stepsize=1):
        return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)

low_t = []
chunk = 30
for day_back in range(-3282, 0, chunk):
    start = DateTime(day_back)
    if (day_back + chunk) < 0:
        stop = DateTime(day_back + chunk)
    else:
        continue
    telem = fetch.get_telem(['AOKALSTR'] + res_msids, start=start, stop=stop,
                            select_events='dwells')
    telem.interpolate(4.1)
    for slot in [1, 2, 5, 6]:
        telem['AORESY{}'.format(slot)] = telem['AORESY{}_1'.format(slot)]
        telem['AORESZ{}'.format(slot)] = telem['AORESZ{}_1'.format(slot)]

    yres_ok =  np.column_stack([np.abs(np.degrees(telem['AORESY{}'.format(slot)].vals)) * 3600 < 20
                                for slot in range(0, 8)])
    zres_ok =  np.column_stack([np.abs(np.degrees(telem['AORESY{}'.format(slot)].vals)) * 3600 < 20
                             for slot in range(0, 8)])
    yres_great =  np.column_stack([np.abs(np.degrees(telem['AORESY{}'.format(slot)].vals)) * 3600 < 5
                             for slot in range(0, 8)])
    zres_great =  np.column_stack([np.abs(np.degrees(telem['AORESY{}'.format(slot)].vals)) * 3600 < 5
                                   for slot in range(0, 8)])

    kalstr = telem['AOKALSTR'].vals.astype(int)
    diff = np.sum((yres_ok != yres_great) | (zres_ok != zres_great), axis=1)
    # find the intervals (with good onboard residual data) where the nonzero
    # diff and low kalstr could have caused a problem.
    for e in consecutive(np.flatnonzero(((kalstr - diff) < 3) & (diff != 0))):
        if len(e) > 8:
            low_t.append(DateTime(telem['AORESY0'].times[e[0]]).date)

