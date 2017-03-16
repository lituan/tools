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


def leaf(labels, similarities, cutoff):

    print 'before leaf ',len(labels)
    matrix = [map(lambda x: 1 if x > cutoff else 0, row)
              for row in similarities]
    for i in range(len(matrix)):
        matrix[i][i] = 0

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
    nr_labels = [i for i in range(len(matrix)) if not i in remove]

    nr_similarities = [similarities[i] for i in range(len(similarities)) if not i in remove]
    nr_similarities = [[row[i] for i in range(
        len(similarities)) if not i in remove] for row in nr_similarities]
    nr_labels = [labels[i] for i in range(len(similarities)) if not i in remove]

    print 'after leaf ',len(nr_labels)
    return nr_labels, nr_similarities


def leaf_seqs(seqs,cutoffs=[0.9]):
    seqnames = [seq[0] for seq in seqs]
    similarities = get_similarity(seqs)
    # back up similarities
    pickle.dump(similarities,open('similarities.pickle','w'))
    nr_seqs_list = []
    for cutoff in cutoffs:
        nr_names,nr_similarities = leaf(seqnames, similarities, cutoff)
        nr_seqs = [seq for seq in seqs if seq[0] in nr_names]
        nr_seqs_list.append((nr_seqs,cutoff))
    return nr_seqs_list

def get_pdb_chain_sequences(p):
    # p format: [(pdbid,chain)...]
    pdbids = ','.join([pdbid for pdbid,chain in p])
    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbids,
        'customReportColumns':'structureId,uniprotAcc,entityId,resolution,chainLength,releaseDate,sequence',
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
    lines = [(line[0],line[1],line[-1]) for line in lines if (line[0],line[1]) in p]
    lines = [line for line in lines if line[-1]]
    return lines


def deredundant(pdb_lig_interactions,cutoffs):
    pdb_lig_chain_interactions = []
    for pdbid,lig,interactions in pdb_lig_interactions:
        chains = []
        for interaction_type,interaction_res in interactions:
            chain = [r.split('_')[0] for r in interaction_res]
            chains += chain
        main_interaction_chain = Counter(chains).most_common()[0][0]
        pdb_lig_chain_interactions.append([pdbid,lig,main_interaction_chain,interactions])

    pdb_chain = [(pdbid,chain) for pdbid,_,chain,_ in pdb_lig_chain_interactions]

    pdb_chain_sequences = get_pdb_chain_sequences(pdb_chain)
    pdb_chain_seqs = [(pdbid+'_'+chain,seq) for pdbid,chain,seq in pdb_chain_sequences]

    pdb_chain_seqs_nr_list = leaf_seqs(pdb_chain_seqs,cutoffs)
    pdb_lig_interactions_nr_list = []
    for pdb_chain_seqs_nr,cutoff in pdb_chain_seqs_nr_list:
        pdb_chain_nr = [pdb_chain for pdb_chain,seq in pdb_chain_seqs_nr]
        pdb_lig_chain_interactions_nr = [p for p in pdb_lig_chain_interactions if p[0]+'_'+p[2] in pdb_chain_nr]

        pdb_lig_interactions_nr = [(p[0],p[1],p[3]) for p in pdb_lig_chain_interactions_nr]
        pdb_lig_interactions_nr_list.append((pdb_lig_interactions_nr,cutoff))

    return pdb_lig_interactions_nr_list


def main():
    cutoff = 0.9
    pdb_lig_interactions = pickle.load(open(sys.argv[-1]))

    fname = '.'.join(os.path.split(sys.argv[-1])[-1].split('.')[:-1])

    cutoffs = [0.4,0.9]

    pdb_lig_interactions_nr_list = deredundant(pdb_lig_interactions,cutoffs)

    for pdb_lig_interactions_nr,cutoff in pdb_lig_interactions_nr_list:

        pickle.dump(pdb_lig_interactions_nr,open(fname+'_nr_'+str(cutoff)+'.pickle','w'))

if __name__ == "__main__":
    main()

