#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
for each lig, there is one main interacting chain
thus, deredundancy is applied to this main interacting chain
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
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist


def align(p):
    s1,s2,seq1,seq2 = p
    if seq1 == seq2:
        return s1,s2,1.0
    matrix = matlist.blosum62
    gap_open = -10  # usual value
    gap_extend = -0.5  # usual value

    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)

    seq1 = alns[0][0]
    seq2 = alns[0][1]
    identity = [1 for i, s in enumerate(seq1) if s == seq2[i]]
    identity = 1.0 * len(identity)/ len(seq1)

    return s1,s2,float('{0:<4.2f}'.format(identity))

def get_similarity(seqs):
    seq_pairs = []
    seq_num = len(seqs)
    for i in range(seq_num):
        for j in range(seq_num):
            if j > i:
                seq_pairs.append((i,j,seqs[i][1],seqs[j][1]))
    p = Pool(6)
    results = p.map(align,seq_pairs)
    p.close()

    results = sorted(results)
    scores = np.ones(shape=(seq_num,seq_num))
    for i,j,s in results:
        scores[i][j] = s
        scores[j][i] = s

    return scores


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


def get_pdb_chain_sequences(p):
    # p format: [(pdbid,chain)...]
    pdbid,inter_id,chain = p
    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbid,
        'customReportColumns':'structureId,sequence',
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

    pdb_chain = [(p[0],p[1],p[-1][0]) for p in pdb_phos_sites]
    p = Pool(6)
    pdb_chain_seqs = p.map(get_pdb_chain_sequences,pdb_chain)
    p.close()

    pdb_chain_seqs = [(pdbid+'_'+inter_id+'_'+chain,seq) for pdbid,inter_id,chain,seq in pdb_chain_seqs]
    similarities = get_similarity(pdb_chain_seqs)
    names  =  [p[0] for p in pdb_chain_seqs]

    nr_list = []
    for cutoff in cutoffs:
        matrix = []
        for i in range(len(pdb_chain_seqs)):
            i_pattern = set(['_'.join(pi[1].split('_')[:3]) for pi in pdb_phos_sites[i][8]])
            i_pattern = sorted([r.split('_')[2] for r in i_pattern])
            matrix_i = []
            for j in range(len(pdb_chain_seqs)):
                j_pattern = set(['_'.join(pi[1].split('_')[:3]) for pi in pdb_phos_sites[j][8]])
                j_pattern = sorted([r.split('_')[2] for r in j_pattern])
                if similarities[i][j] >= cutoff and i_pattern == j_pattern:
                    matrix_i.append(1)
                else:
                    matrix_i.append(0)
            matrix.append(matrix_i)

        for i in range(len(matrix)):
            matrix[i][i] = 0

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

