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

MP_DIR = '/data/mpcrit1/mplogs'


def parse_clgps(file):

    # assume this isn't a replan
    summary = {'replan': False,
               'replan_cmds': '',
               'continuity_cmds': '',
               'bcf_cmd_count': 0}

    rawloads = []

    # split the file into paragraphs based on the asterisk separator lines
    chunks = open(file).read().split("\n************")

    # do the yucky, but maintainable, regex parsing of each paragraph
    for piece in chunks:
        if re.search("CONTINUITY/REPLAN", piece):
            if re.search("Replan/Reopen", piece):
                summary['replan'] = True
                run_dir = re.search("Run\sDirectory:\s\S+\/(C\d{3}:\d{4})", piece)
                if run_dir:
                    summary['replan_cmds'] = run_dir.group(1)
                cont_dir = re.search("Continuity\sDirectory:\s\S+\/(C\d{3}:\d{4})", piece)
                if cont_dir:
                    summary['continuity_cmds'] = cont_dir.group(1)
            else:
                cont_dir = re.search("Continuity\sDirectory:\s\S+\/(C\d{3}:\d{4})", piece)
                if cont_dir:
                    summary['continuity_cmds'] = cont_dir.group(1)
        if re.search("\*+PROCESSING\sVALUES\*+", piece):
            exec_time = re.search("EXECUTION\sBEGIN\sTIME:(\S+)", piece)
            if exec_time:
                summary['execution_tstart'] = exec_time.group(1)
            proc_start = re.search("START\sTIME:\s+(\S+)", piece)
            if proc_start:
                summary['processing_tstart'] = proc_start.group(1)
            proc_stop = re.search("STOP\sTIME:\s+(\S+)", piece)
            if proc_stop:
                summary['processing_tstop'] = proc_stop.group(1)
        if re.search("\*+\sINPUT\sREPLAN\sBCF\sINFORMATION\s\*+", piece):
            bcf = re.search("TOTAL\sNUMBER\sOF\sCOMMANDS\sREAD\s=\s(\d+)", piece)
            if bcf:
                summary['bcf_cmd_count'] = int(bcf.group(1))
        if re.search("\*+\sACTUAL\sPLANNING\sPERIOD\sSTART\/STOP\sTIMES\s\*+", piece):
            actual_start = re.search("ACTUAL\sSTART\sTIME:\s(\S+)", piece)
            if actual_start:
                summary['planning_tstart'] = actual_start.group(1)
                summary['year'] = int(summary['planning_tstart'][0:4])
            actual_stop = re.search("ACTUAL\sSTOP\s\sTIME:\s(\S+)", piece)
            if actual_stop:
                summary['planning_tstop'] = actual_stop.group(1)
        if re.search("\*+LOAD\sGENERATED\*+", piece):
            load = {}
            load_seg = re.search("Load\sname:\s+(CL\d{3}:\d{4})", piece)
            if load_seg:
                load['load_segment'] = load_seg.group(1)
            scs = re.search("SCS\sNumber:\s+(\d{3}(\s*\/\s*\d{3})?)", piece)
            if scs:
                load['load_scs'] = scs.group(1)
            first_cmd = re.search("First\sCommand\sTime:\s+(\S+)", piece)
            if first_cmd:
                load['first_cmd_time'] = first_cmd.group(1)
                load['year'] = int(load['first_cmd_time'][0:4])
            last_cmd = re.search("Last\sCommand\sTime:\s+(\S+)", piece)
            if last_cmd:
                load['last_cmd_time'] = last_cmd.group(1)
            rawloads.append(load)

    loads = []
    for rawload in rawloads:
        scs_pat = re.match("^(\d{3})\s*\/\s*(\d{3})$", rawload['load_scs'])
        if scs_pat:
            vehicle_scs = scs_pat.group(1)
            observing_scs = scs_pat.group(2)
            veh_load = rawload.copy()
            veh_load['load_scs'] = vehicle_scs
            obs_load = rawload.copy()
            obs_load['load_scs'] = observing_scs
            loads.append(veh_load)
            loads.append(obs_load)
        else:
            loads.append(rawload)

    return summary, loads


def get_loads():
    week_glob = '[A-Z][A-Z][A-Z][0-9][0-9][0-9][0-9]';
    sum_glob = 'C[0-9][0-9][0-9]?[0-9][0-9][0-9][0-9].sum';
    now = DateTime()

    files = []
    for year in range(2008, int(now.frac_year) + 1):
        cmd = ("find -L {mp_dir}/{year:04d}/{week_glob}/ofls?/mps/ ".format(mp_dir=MP_DIR, year=year, week_glob=week_glob)
               + ' -maxdepth 1 -wholename \"*/ofls?/mps/{sum_glob}\" '.format(sum_glob=sum_glob))
        files.extend(bash(cmd))

    summaries = []
    loads = []
    for filename in files:
        name = re.match(MP_DIR + "/(((\d{4})/(\w{3}\d{4})/ofls(\w))/mps/(C.*\.sum))",
                        filename)
        if name is None:
            continue
        if name and name.groups(4) == 'x' or name.groups(4) == 't':
            continue
        summary, week_loads = parse_clgps(filename)
        summary['file'] = name.group(1)
        summary['shortfile'] = name.group(6)
        summary['mp_dir'] = name.group(2)
        for load in week_loads:
            load['file'] = name.group(1)
            load['mp_dir'] = name.group(2)
            loads.append(load)
        summaries.append(summary)

    loads = Table(loads)
    loads.sort('first_cmd_time')
    summaries = Table(summaries)
    summaries.sort('processing_tstart')
    return summaries, loads


def get_load_segments():
    with Ska.DBI.DBI(dbi='sybase', server='sybase', database='aca', user='aca_read') as db:
        load_segments = db.fetchall("select * from load_segments")
    load_segments = Table(load_segments)
    return load_segments


def get_catalog(obsid, obstime):
    idx = np.searchsorted(load_segments['datestart'], obstime) - 1
    # if the start of the dwell isn't captured in this load segment, skip the dwell
    if not ((load_segments[idx]['datestart'] < obstime)
            & (load_segments[idx]['datestop'] > obstime)):
        print "Skipping {} {}".format(obstime, obsid)
    if load_segments[idx]['load_segment'] == 'CL344:1903_1':
        load_segments[idx]['load_segment'] = 'CL344:1903'
    if load_segments[idx]['load_segment'] == 'CL208:1800_1':
        load_segments[idx]['load_segment'] = 'CL208:1800'
    if load_segments[idx]['load_segment'] == 'CL344:2303_2':
        load_segments[idx]['load_segment'] = 'CL344:2303'
    load_match = ((loads['year'] == load_segments[idx]['year'])
                  & (loads['load_segment'] == load_segments[idx]['load_segment']))

    load = loads[load_match][0]
    summary = summaries[summaries['file'] == load['file']][0]
    mp_dir = summary['mp_dir']
    if summary['bcf_cmd_count'] != 0:
        if (summary['processing_tstart'] > obstime):
            re_match = re.match("C(\d{3}):(\d{4})", summary['replan_cmds'])
            if re_match:
                sum_match = summaries[summaries['shortfile'] == "C{}_{}.sum".format(
                        re_match.group(1),
                        re_match.group(2))]
                if len(sum_match):
                    mp_dir = sum_match[0]['mp_dir']
                else:
                    raise ValueEarror
    cat = get_starcheck_catalog(obsid, mp_dir="/{}/".format(mp_dir))
    if not len(cat['cat']):
        raise ValueError
    targdelta = DateTime(obstime).secs - DateTime(cat['manvr'][-1]['mp_targquat_time']).secs
    return cat


def get_pcad(dwell):
    msids = ['AOACASEQ', 'AOPCADMD', 'AOKALSTR', 'AONSTARS',
             'AOATTQT1', 'AOATTQT2', 'AOATTQT3', 'AOATTQT4', 'AOACIMSS']
    slots = range(0, 8)
    slot_cols = ['AOACMAG', 'AOACYAN', 'AOACZAN', 'AOACFCT', 'AOIMAGE',
                 'AOACIDP', 'AOACIIR', 'AOACIMS', 'AOACISP']
    slot_msids = ["{}{}".format(field, slot) for field in slot_cols
                  for slot in slots]
    msids.extend(slot_msids)
    return fetch.MSIDset(msids, dwell.start, dwell.stop)


if 'loads' not in globals():
    summaries, loads = get_loads()
if 'load_segments' not in globals():
    load_segments = get_load_segments()


import agasc
from Ska.quatutil import radec2eci

orig_bad_obsids = []
bad_obsids = []

def consecutive(data, stepsize=1):
        return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)

ds = events.dwells.filter(start='2008:007')
#ds = events.dwells.filter(obsid=17321)
low_obsids = []
for d in ds:
    obsid = d.get_obsid()
    print "obsid {} start {}".format(obsid, d.manvr.start)
    pcad_data = get_pcad(d)
    q_atts = Quat(normalize(np.column_stack([pcad_data['AOATTQT1'].vals,
                                             pcad_data['AOATTQT2'].vals,
                                             pcad_data['AOATTQT3'].vals,
                                             pcad_data['AOATTQT4'].vals])))
    try:
        cat = get_catalog(obsid, d.manvr.start)
    except (ValueError, IndexError) as e:
        print "Skipping {} {}".format(obsid, e)
        continue
    gcat = cat['cat'][(cat['cat']['type'] == 'BOT') | (cat['cat']['type'] == 'GUI')]
    yag_offs = np.zeros((len(pcad_data['AOKALSTR'].vals), len(gcat)))
    zag_offs = np.zeros((len(pcad_data['AOKALSTR'].vals), len(gcat)))
    slot_ok_old = {}
    slot_ok_new = {}
    for idx, entry in enumerate(gcat):
        slot = entry['slot']
        #ok = pcad_data['AOACFCT{}'.format(slot)].vals == 'TRAK'
        star = agasc.get_star(entry['id'], date=d.manvr.start)
        eci = radec2eci(star['RA_PMCORR'], star['DEC_PMCORR'])
        d_aca = np.dot(q_atts.transform.transpose(0,2,1), eci)
        yag = np.degrees(np.arctan2(d_aca[:, 1], d_aca[:, 0])) * 3600
        zag = np.degrees(np.arctan2(d_aca[:, 2], d_aca[:, 0])) * 3600
        yag_offs[:, idx] = yag - pcad_data['AOACYAN{}'.format(entry['slot'])].vals
        zag_offs[:, idx] = zag - pcad_data['AOACZAN{}'.format(entry['slot'])].vals
    diff = np.sum(((np.abs(yag_offs) > 5) & (np.abs(yag_offs) < 20))
                  | ((np.abs(zag_offs) > 5) & (np.abs(zag_offs) < 20)), axis=1)
    kalstr = pcad_data['AOKALSTR'].vals.astype(int)
    for e in consecutive(np.flatnonzero(((kalstr - diff) < 3) & (diff != 0))):
        if len(e) > 8:
            low_obsids.append(obsid)

low_obsids = np.unique(low_obsids)


#        rough_kalstr = np.sum([slot_ok_new[i].astype(int) for i in guide_stars], axis=0)
#        def consecutive(data, stepsize=1):
#            return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)
#        extents = consecutive(np.flatnonzero(rough_kalstr < 2))
#        for extent in extents:
#            if len(extent) > 30:
#                print "kalstr error on {}".format(obsid)
#                bad_obsids.append({'obsid': obsid,
#                                   'time': pcad_data['AOKALSTR'].times[extent[0]]})
#        orig_extents = consecutive(np.flatnonzero(pcad_data['AOKALSTR'].vals.astype(int) < 2))
#        for extent in orig_extents:
#            if len(extent) > 30:
#                print "orig kalstr error on {}".format(obsid)
#                orig_bad_obsids.append({'obsid': obsid,
#                                        'time': pcad_data['AOKALSTR'].times[extent[0]]})
##    raise ValueError
#        except:
#            continue
#    if len(bad_obsids):
#        Table(bad_obsids).write("bad_obsids.dat", format='ascii')
#    if len(orig_bad_obsids):
#        Table(orig_bad_obsids).write("orig_bad_obsids.dat", format='ascii')
