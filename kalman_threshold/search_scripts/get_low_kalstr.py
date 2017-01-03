from kadi import events
from mica.quaternion import Quat, normalize
from mica.starcheck import get_starcheck_catalog
import os
import time
import re
from Chandra.Time import DateTime
import Ska.DBI
from Ska.Shell import bash
from astropy.table import Table
import numpy as np
from Ska.engarchive import fetch

kalstr_obsids = []
for year in range(2008, 2017):
    dat = fetch.Msid('AOKALSTR', '{}:001'.format(year), '{}:001'.format(year+1))
    lowkals = dat.logical_intervals('<=', '2 ')
    lowkals = lowkals[lowkals['duration'] < 120]  # too-long intervals are spurious
    bad = lowkals['duration'] > 60
    lowkals = lowkals[bad]
    for interval in lowkals:
        ds = events.dwells.filter(start=interval['tstart'], stop=interval['tstop'])
        kalstr_obsids.extend([d.get_obsid() for d in ds])

kalstr_obsids = np.unique(kalstr_obsids)
kalstr_obsids = kalstr_obsids[kalstr_obsids != 0]
