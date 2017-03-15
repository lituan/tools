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


def leaf(labels, matrix):

    print 'before leaf ',len(labels)

    adjlist = [[i for i,n in enumerate(row ) if n] for row in matrix]
    neighbors = []
    remove = []
    neighbors = [set(n) for i, n in enumerate(adjlist)]
    # detect possible max clique
    max_neighbors = max(len(l) for l in neighbors)
    # the possible clique size is 2 to max_neighbors+1, so the possible
    # neighborsize is 1 to max_neighbors
    for clique_num in range(1, max_neighbors + 1):
        nodes_index = set([i for i, l in enumerate(
            neighbors) if len(l) == clique_num])
        for i in nodes_index:
            if not i in remove:  # do not compute removed vertex
                # a clique is set of vertex connecting to each other
                nodesofinterest = neighbors[i].union([i])
                if set.intersection(*[neighbors[i].union([i]) for i in nodesofinterest]) == nodesofinterest:
                    # detect vertex without linking to outside vertex
                    in_clique = [i for i in nodesofinterest if not neighbors[
                        i].union([i]).difference(nodesofinterest)]
                    # keep one of the vertex without linking to outside vertex,
                    # remove rest
                    if in_clique:
                        keep = [in_clique[0]]
                        remove_iter = nodesofinterest.difference(set(keep))
                        for r in remove_iter:
                            if not r in remove:  # do not compute removed vertex
                                for i in range(len(neighbors)):
                                    if r in neighbors[i]:
                                        neighbors[i].remove(r)
                        remove += remove_iter

    nr_matrix = [matrix[i] for i in range(len(matrix)) if not i in remove]
    nr_matrix = [[row[i] for i in range(
        len(matrix)) if not i in remove] for row in nr_matrix]
    nr_labels = [i for i in range(len(matrix)) if not i in remove]
    # continue to remove the one with most neighbors until no vertex has
    # neighbors, removed vertex is not considered
    while max([len(r) for i, r in enumerate(neighbors) if not i in remove]) > 0:

        max_index = max([(len(r), i) for i, r in enumerate(neighbors) if not i in remove])[1]
        remove.append(max_index)
        for i in set(range(len(neighbors))).difference(set(remove)): # do not compute remove vertex
            if max_index in neighbors[i]:
                neighbors[i].remove(max_index)

    nr_matrix = [matrix[i] for i in range(len(matrix)) if not i in remove]
    nr_matrix = [[row[i] for i in range(
        len(matrix)) if not i in remove] for row in nr_matrix]
    nr_labels = [labels[i] for i in range(len(matrix)) if not i in remove]

    nr_clusters = [(i,a) for i,a in enumerate(adjlist) if not i in remove]
    nr_labels_clusters = [([labels[i]]+[labels[ai] for ai in a]) for i,a in nr_clusters]

    print 'after leaf ',len(nr_labels)
    return nr_labels, nr_labels_clusters


def get_pdb_chain_clusternum(p):
    # p format: [(pdbid,chain)...]
    pdbid,inter_id,chain,cutoff = p
    if cutoff == 1.0:
        clusternum = 'clusterNumber100'
    elif cutoff == 0.95:
        clusternum = 'clusterNumber95'
    elif cutoff == 0.90:
        clusternum = 'clusterNumber90'
    elif cutoff == 0.70:
        clusternum = 'clusterNumber70'
    elif cutoff == 0.50:
        clusternum = 'clusterNumber50'
    elif cutoff == 0.40:
        clusternum = 'clusterNumber40'
    elif cutoff == 0.30:
        clusternum = 'clusterNumber30'
    else:
        print 'wrong cutoff'
        print 'permitted cutoffs: 0.3,0.4,0.5,0.7,0.9,0.95'
        sys.exit()
    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbid,
        'customReportColumns':'structureId,'+clusternum,
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
    lines = [(line[0],line[1],line[-1]) for line in lines if line[1] == chain][0]
    return lines[0],inter_id,lines[1],lines[-1]


def deredundant(pdb_phos_sites,cutoffs):

    nr_list = []
    for cutoff in cutoffs:
        pdb_chain = [(p[0],p[1],p[-1][0],cutoff) for p in pdb_phos_sites]
        p = Pool(6)
        pdb_chain_clusternums = p.map(get_pdb_chain_clusternum,pdb_chain)
        p.close()

        result = []
        get_pdb_chain_clusternum(pdb_chain[0])
        sys.exit()

        clusternums = {}
        for pdbid,inter_id,chain,clusternum in pdb_chain_clusternums:
            clusternums[pdbid+'_'+inter_id+'_'+chain] = clusternum

        matrix = []
        for i in range(len(pdb_chain)):
            i_pattern = set(['_'.join(pi[1].split('_')[:3]) for pi in pdb_phos_sites[i][8]])
            i_pattern = sorted([r.split('_')[2] for r in i_pattern])
            ic =  pdb_phos_sites[i][0] + '_' + pdb_phos_sites[i][1] + '_' + pdb_phos_sites[i][-1][0]
            matrix_i = []
            for j in range(len(pdb_chain)):
                j_pattern = set(['_'.join(pi[1].split('_')[:3]) for pi in pdb_phos_sites[j][8]])
                j_pattern = sorted([r.split('_')[2] for r in j_pattern])
                jc =  pdb_phos_sites[j][0] + '_' + pdb_phos_sites[j][1] + '_' + pdb_phos_sites[j][-1][0]
                if clusternums[ic] == clusternums[jc] and i_pattern == j_pattern:
                    matrix_i.append(1)
                else:
                    matrix_i.append(0)
            matrix.append(matrix_i)

        for i in range(len(matrix)):
            matrix[i][i] = 0

        names  =  [p[0]+'_'+p[1]+'_'+p[2] for p in pdb_chain]
        nr_names,nr_names_clusters = leaf(names,matrix)
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

