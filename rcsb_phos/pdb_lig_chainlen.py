#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
check chain length of pdb chains
"""

import os
import sys
import urllib
import urllib2
import cPickle as pickle
import numpy as np
from collections import OrderedDict
from collections import Counter
from multiprocessing import Pool

def get_pdb_chain_length(p):
    # p format: [(pdbid,chain)...]
    pdbid,lig,lig_chain,interacting_chains = p
    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbid,
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
    lig_chain_len = [line[5] for line in lines if line[1] == lig_chain]
    if lig_chain_len:
        lig_chain_len = lig_chain_len[0]
    interacting_chains_len = [(line[1],line[5]) for line in lines if line[1] in interacting_chains]
    if lig_chain_len and interacting_chains_len:
        return (pdbid,lig,lig_chain_len,interacting_chains_len)

def main():
    pdb_lig_interactions = pickle.load(open(sys.argv[-1]))

    pdb_lig_chain = []
    for pdbid,lig,interactions in pdb_lig_interactions:
        lig_chain = lig.split('_')[0]
        interacting_chains = []
        for interaction_type,interaction_res in interactions:
            chain = [r.split('_')[0] for r in interaction_res]
            interacting_chains += chain
        interacting_chains = set(interacting_chains)
        pdb_lig_chain.append((pdbid,lig,lig_chain,interacting_chains))

    result = []
    for p in pdb_lig_chain:
        result.append(get_pdb_chain_length(p))

    fname = '.'.join(os.path.split(sys.argv[-1])[-1].split('.')[:-1])
    pickle.dump(result,open(fname+'_chain_len.pickle','w'))
    with open(fname+'_chain_len.txt','w') as w_f:
        for pdb,lig,lig_chain_len,interacting_chains_len in result:
            print >> w_f,'{0:<10}{1:<15}{2:<10}{3:<}'.format(pdb,lig,lig_chain_len,interacting_chains_len)

if __name__ == "__main__":
    main()

