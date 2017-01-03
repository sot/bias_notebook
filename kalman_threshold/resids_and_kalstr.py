import numpy as np

from mica.quaternion import Quat, normalize
from Ska.quatutil import radec2eci
from Chandra.Time import DateTime
import agasc
from Ska.engarchive import fetch


def get_telem(dwell):
    subformat = fetch.Msid('COTLRDSF', dwell.start, dwell.stop)
    msids = ['AOACASEQ', 'AOPCADMD', 'AOKALSTR',  'COTLRDSF', 'COBSRQID',
             'AOATTQT1', 'AOATTQT2', 'AOATTQT3', 'AOATTQT4', 'CORADMEN', '3TSCPOS', '3TSCMOVE']
    # don't bother getting the residuals if the dwell has PCAD subformat
    if not np.any(subformat.vals == 'PCAD'):
        res_msids = ['AORESY0', 'AORESY1_1', 'AORESY2_1', 'AORESY3', 'AORESY4',
                     'AORESY5_1', 'AORESY6_1', 'AORESY7',
                     'AORESZ0', 'AORESZ1_1',
                     'AORESZ2_1', 'AORESZ3', 'AORESZ4', 'AORESZ5_1', 'AORESZ6_1',
                     'AORESZ7']
        msids = msids + res_msids
    slot_cols = ['AOACFID', 'AOACMAG', 'AOACYAN', 'AOACZAN', 'AOACFCT', 'AOIMAGE',
                 'AOACIDP', 'AOACIIR', 'AOACIMS', 'AOACISP']
    slots = range(0, 8)
    slot_msids = ["{}{}".format(field, slot) for field in slot_cols
                  for slot in slots]
    msids = msids + slot_msids
    if dwell.start > '2015:251':
        pcad_data = fetch.Msidset(msids + ['AOACIMSS'], start=dwell.start, stop=dwell.stop)
    else:
        pcad_data = fetch.Msidset(msids, start=dwell.start, stop=dwell.stop)

    if 'AORESY0' in pcad_data:
        for slot in [1, 2, 5, 6]:
            pcad_data['AORESY{}'.format(slot)] = pcad_data['AORESY{}_1'.format(slot)]
            pcad_data['AORESZ{}'.format(slot)] = pcad_data['AORESZ{}_1'.format(slot)]

    return pcad_data


def kal(dwell, telem, limit=20, catalog=None, nowflags=False):

    cat = catalog

    # Track status
    fids = np.column_stack([(telem['AOACFID{}'.format(slot)].vals == 'FID ')
                            for slot in range(0, 8)])
    trak = np.column_stack([(telem['AOACFCT{}'.format(slot)].vals == 'TRAK')
                            for slot in range(0, 8)])
    # Flags
    ir = np.column_stack([(telem['AOACIIR{}'.format(slot)].vals == 'OK ')
                          for slot in range(0, 8)])
    sp = np.column_stack([(telem['AOACISP{}'.format(slot)].vals == 'OK ')
                          for slot in range(0, 8)])
    dp_date = DateTime('2013:297:11:25:52.000').secs
    dp = np.column_stack([((telem['AOACIDP{}'.format(slot)].vals == 'OK ')
                           | (telem['AOKALSTR'].times > dp_date))
                          for slot in range(0, 8)])
    if dwell.start > '2015:251':
        # I'm not sure about the fetch grid if we use fetch interpolate, so just use
        # a sorted search to see if the MSS flag should apply
        mss = telem['AOACIMSS'].vals == 'ENAB'
        mss_times = telem['AOACIMSS'].times
        mss_at_times = mss[np.searchsorted(mss_times[1:-1], telem['AOACIMS0'].times) - 1]
        ms = np.column_stack([((telem['AOACIMS{}'.format(slot)].vals == 'OK ')
                               | ~mss_at_times)
                              for slot in range(0, 8)])
    else:
        ms = np.column_stack([(telem['AOACIMS{}'.format(slot)].vals == 'OK ')
                              for slot in range(0, 8)])

    # use rolled-by-4 for ~last 4.1 sample
    last_trak = np.roll(trak, 4, axis=0)
    last_trak[0] = True

    # Calc centroid residuals using CYAN/CZAN
    q_atts = Quat(normalize(np.column_stack([telem['AOATTQT1'].vals,
                                             telem['AOATTQT2'].vals,
                                             telem['AOATTQT3'].vals,
                                             telem['AOATTQT4'].vals])))

    # guide and bot slots
    gcat = cat[(cat['type'] == 'BOT') | (cat['type'] == 'GUI')]
    # make a couple of structures for the offsets
    yag_offs = np.zeros((len(telem['AOKALSTR'].vals), 8))
    zag_offs = np.zeros((len(telem['AOKALSTR'].vals), 8))
    for idx, entry in enumerate(gcat):
        slot = entry['slot']
        #ok = telem['AOACFCT{}'.format(slot)].vals == 'TRAK'
        star = agasc.get_star(entry['id'], date=dwell.manvr.start)
        eci = radec2eci(star['RA_PMCORR'], star['DEC_PMCORR'])
        d_aca = np.dot(q_atts.transform.transpose(0,2,1), eci)
        yag = np.degrees(np.arctan2(d_aca[:, 1], d_aca[:, 0])) * 3600
        zag = np.degrees(np.arctan2(d_aca[:, 2], d_aca[:, 0])) * 3600
        yag_offs[:, slot] = yag - telem['AOACYAN{}'.format(entry['slot'])].vals
        zag_offs[:, slot] = zag - telem['AOACZAN{}'.format(entry['slot'])].vals

    if nowflags:
        kal = (~fids & trak & last_trak & ir & sp
                & (np.abs(yag_offs) < limit) & (np.abs(zag_offs) < limit))
    else:
        kal = (~fids & trak & last_trak & ir & sp & dp & ms
                & (np.abs(yag_offs) < limit) & (np.abs(zag_offs) < limit))

    raise ValueError
    return telem['AOKALSTR'].times, kal
