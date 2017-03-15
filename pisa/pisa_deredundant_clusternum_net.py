#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
for each lig, there is one main interacting chain
thus, deredundancy is applied to this main interacting chain

use pdb clusternum to detect similarity
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


def get_pdb_chain_clusternum(p):
    # p format: [(pdbid,chain)...]
    pdbid,inter_id,chain,cutoff = p
    if cutoff == 1.0:
        clusternum = 'clusterNumber100'
        ranknum = 'rankNumber100'
    elif cutoff == 0.95:
        clusternum = 'clusterNumber95'
        ranknum = 'rankNumber95'
    elif cutoff == 0.90:
        clusternum = 'clusterNumber90'
        ranknum = 'rankNumber90'
    elif cutoff == 0.70:
        clusternum = 'clusterNumber70'
        ranknum = 'rankNumber70'
    elif cutoff == 0.50:
        clusternum = 'clusterNumber50'
        ranknum = 'rankNumber50'
    elif cutoff == 0.40:
        clusternum = 'clusterNumber40'
        ranknum = 'rankNumber40'
    elif cutoff == 0.30:
        clusternum = 'clusterNumber30'
        ranknum = 'rankNumber30'
    else:
        print 'wrong cutoff'
        print 'permitted cutoffs: 0.3,0.4,0.5,0.7,0.9,0.95'
        sys.exit()
    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbid,
        'customReportColumns':'structureId,'+clusternum+','+ranknum,
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
    lines = [(line[0],line[1],int(line[2]),int(line[3])) for line in lines if line[1] == chain][0]
    return lines[0],inter_id,lines[1],lines[2],lines[3]


def deredundant(pdb_phos_sites,cutoffs):

    nr_list = []
    for cutoff in cutoffs:
        pdb_chain = [(p[0],p[1],p[-1][0],cutoff) for p in pdb_phos_sites]
        p = Pool(6)
        pdb_chain_clusternums = p.map(get_pdb_chain_clusternum,pdb_chain)
        p.close()

        clusternums = OrderedDict()
        for i,piccr in enumerate(pdb_chain_clusternums):
            pdbid,inter_id,chain,clusternum,ranknum = piccr
            ic = pdbid+'_'+inter_id+'_'+chain
            i_pattern = set(['_'.join(pi[1].split('_')[:3]) for pi in pdb_phos_sites[i][8]])
            i_pattern = sorted([r.split('_')[2] for r in i_pattern])
            i_pattern = '_'.join(i_pattern)
            cc = str(clusternum ) + '-' + i_pattern
            if not cc in clusternums.keys():
                clusternums[cc] = [(ranknum,ic)]
            else:
                clusternums[cc].append((ranknum,ic))

        nr_names = [sorted(v,key=lambda x: x[0],reverse=False)[0][1] for k,v in clusternums.iteritems() ]
        nr_names_clusters = [[vi[1] for vi in sorted(v,key=lambda x: x[0],reverse=False)] for k,v in clusternums.iteritems()]

        print cutoff,len(nr_names)

        pdb_phos_sites_nr = [p for p in pdb_phos_sites if p[0]+'_'+p[1]+'_'+p[-1][0] in nr_names]
        pdb_phos_sites_nr_clusters= [[p for p in pdb_phos_sites if p[0]+'_'+p[1]+'_'+p[-1][0] in cluster ] for cluster in nr_names_clusters]
        nr_list.append([pdb_phos_sites_nr,pdb_phos_sites_nr_clusters,cutoff])

    return nr_list


def main():
    pdb_phos_sites = pickle.load(open(sys.argv[-1]))

    fname = '.'.join(os.path.split(sys.argv[-1])[-1].split('.')[:-1])

    cutoffs = [0.3,0.9,1.0]

    phos_sites_nr_list= deredundant(pdb_phos_sites,cutoffs)

    for phos_sites,clusters,cutoff in phos_sites_nr_list:

        pickle.dump(phos_sites,open(fname+'_nr_'+str(cutoff)+'.pickle','w'))

        with open(fname+'_nr_'+str(cutoff)+'_cluster.txt','w') as w_f:
            for c in clusters:
                for ci in c:
                    print >> w_f,ci
                print >> w_f,''


if __name__ == "__main__":
    main()

