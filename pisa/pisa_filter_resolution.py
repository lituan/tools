#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import urllib
import urllib2
import cPickle as pickle
from multiprocessing import Pool

"""
filter non_bond, low-resolution, same-chain, non-phos, two-phos interfaces
"""
def write_pdb_interfaces(pdb_interfaces,fname):
    with open(fname+'.txt','w') as w_f:
        for i in pdb_interfaces:
            pdbid,interface_id,interface_area,inter_solvent,pvalue,binding_energy,interaction_res,chains = i
            print >> w_f, '{0:<8}{1:<8}{2:<8.2f}{3:<8.2f}{4:<8.2f}{5:<8.2f}{6:<}\t{7:<}'.format(pdbid,interface_id,interface_area,inter_solvent,pvalue,binding_energy,interaction_res,chains)
        pickle.dump(pdb_interfaces,open(fname+'.pickle','w'))


def filter_resolution(pdb_interfaces,cutoff):

    pdbids = [p[0] for p in pdb_interfaces]
    pdbids = ','.join(pdbids)
    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbids,
        'customReportColumns':'structureId,uniprotAcc,entityId,resolution,chainLength,releaseDate',
        'service':'wsfile',
        'format':'csv',
    }
    data = urllib.urlencode(data)
    req = urllib2.Request(url,data)
    response = urllib2.urlopen(req)
    lines = response.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line for line in lines if line]
    lines = [line.split(',') for line in lines]
    lines = [[w.strip('"') for w in line] for line in lines]

    # filter pdbs with lower resolution
    lines = [line for line in lines[1:] if line[4]]
    pdbids = [line[0] for line in lines if float(line[4]) <= cutoff]

    good_pdb_interfaces = [p for p in pdb_interfaces if p[0] in pdbids]
    return good_pdb_interfaces


def main():
    pdb_interfaces = pickle.load(open(sys.argv[-1]))


    cutoffs = [1.0,2.0,2.5,3.0,3.5]
    for resolution in cutoffs:
        good_pdb_interfaces = filter_resolution(pdb_interfaces,resolution)
        print resolution, len(good_pdb_interfaces)
        write_pdb_interfaces(good_pdb_interfaces,'filter_low_resolution_'+str(resolution))



if __name__ == "__main__":
    main()

