#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
use leaf alogrithm to remove redundancy
this algorithm is published in this following paper
1. Bull, S. C., Muldoon, M. R. & Doig, A. J. Maximising the Size of Non-Redundant Protein Datasets Using Graph Theory. PLoS One 8, (2013).
Simon Bull gives an implementation for python 3, hosted at https://github.com/SimonCB765/Leaf
this implementation is for python 2.7

this script can read in a multiple alignment file in fasta format, compute pairwise similarities, and remove redundancy accordint to similarity cutoff
you can choose to use igraph to plot the network
usage python leaf.py test.fa
"""
import sys
import os
import igraph


def read_msa(msa_f):
    with open(msa_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines if line]
        pro_line_num = [i for i, line in enumerate(
            lines) if '>' in line] + [len(lines)]
        seqs = [lines[n:pro_line_num[i + 1]]
                for i, n in enumerate(pro_line_num[:-1])]
        seqs = [(seq[0].split()[0][1:], ''.join(seq[1:])) for seq in seqs]

        return seqs


def get_pim(seqs):

    def pim(seq1, seq2):
        identity = len([i for i, s in enumerate(seq1) if s == seq2[i]])
        return identity * 1.0 / len(seq1)

    scores = []
    seqlen = len(seqs)
    for i in range(seqlen):
        score_i = []
        for j in range(seqlen):
            if j < i:
                score_i.append(scores[j][i])
            elif j > i:
                score_i.append(pim(seqs[i][1], seqs[j][1]))
            else:
                score_i.append(1.0)
        scores.append(score_i)
    return scores

def leaf(labels, similarities, cutoff, filename):

    matrix = [map(lambda x: 1 if x > cutoff else 0, row)
              for row in similarities]
    for i in range(len(matrix)):
        matrix[i][i] = 0

    # use igraph to plot the initial network
    graph = igraph.Graph.Adjacency(matrix, mode='undirected')
    igraph.plot(graph, filename + '.png', vertex_label=range(len(labels)))

    adjlist = [[i for i,n in enumerate(row ) if n] for row in matrix]
    neighbors = []
    remove = []
    for i,a in enumerate(adjlist):
        print '{0}:{1},'.format(i,a)
    # transform adjlist to set
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
                print 'initial nodesofinterest: ',nodesofinterest
                if set.intersection(*[neighbors[i].union([i]) for i in nodesofinterest]) == nodesofinterest:
                    print 'clique nodesofinterest: ',nodesofinterest
                    # detect vertex without linking to outside vertex
                    in_clique = [i for i in nodesofinterest if not neighbors[
                        i].union([i]).difference(nodesofinterest)]
                    # keep one of the vertex without linking to outside vertex,
                    # remove rest
                    if in_clique:
                        print 'in_clique: ',in_clique
                        keep = [in_clique[0]]
                        print 'keep: ',keep
                        remove_iter = nodesofinterest.difference(set(keep))
                        print 'remove_iter: ',remove_iter
                        for r in remove_iter:
                            if not r in remove:  # do not compute removed vertex
                                print 'remove: ',r
                                for i in range(len(neighbors)):
                                    if r in neighbors[i]:
                                        neighbors[i].remove(r)
                        remove += remove_iter

    print 'after leaf: ',neighbors

    nr_matrix = [matrix[i] for i in range(len(matrix)) if not i in remove]
    nr_matrix = [[row[i] for i in range(
        len(matrix)) if not i in remove] for row in nr_matrix]
    graph = igraph.Graph.Adjacency(nr_matrix, mode='undirected')
    nr_labels = [i for i in range(len(matrix)) if not i in remove]
    igraph.plot(graph, filename + '_leaf.png', vertex_label=nr_labels)
    # continue to remove the one with most neighbors until no vertex has
    # neighbors, removed vertex is not considered
    while max([len(r) for i, r in enumerate(neighbors) if not i in remove]) > 0:

        max_index = max([(len(r), i) for i, r in enumerate(neighbors) if not i in remove])[1]
        print 'remove: ',max_index
        remove.append(max_index)
        for i in set(range(len(neighbors))).difference(set(remove)): # do not compute remove vertex
            if max_index in neighbors[i]:
                neighbors[i].remove(max_index)

    print 'final remove: ',remove

    nr_matrix = [matrix[i] for i in range(len(matrix)) if not i in remove]
    nr_matrix = [[row[i] for i in range(
        len(matrix)) if not i in remove] for row in nr_matrix]
    nr_labels = [i for i in range(len(matrix)) if not i in remove]

    # plot non-redundant notwork
    graph = igraph.Graph.Adjacency(nr_matrix, mode='undirected')
    igraph.plot(graph, filename + '_nr.png', vertex_label=nr_labels)

    nr_similarities = [similarities[i] for i in range(len(similarities)) if not i in remove]
    nr_similarities = [[row[i] for i in range(
        len(similarities)) if not i in remove] for row in nr_similarities]
    nr_labels = [labels[i] for i in range(len(similarities)) if not i in remove]

    return nr_labels, nr_similarities


def main():
    seqs = read_msa(sys.argv[-1])
    filename = os.path.splitext(os.path.split(sys.argv[-1])[1])[0]
    seqnames = [seq[0] for seq in seqs]
    similarities = get_pim(seqs)
    for cutoff in [0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]:
    # for cutoff in [0.8]:
        nr_names,nr_similarities = leaf(seqnames, similarities, cutoff, filename+'_'+str(cutoff))
        nr_seqs = [seq for seq in seqs if seq[0] in nr_names]
        with open(filename+'_nr_seqs_'+str(cutoff)+'.fas','w') as w_f:
            for pro,seq in nr_seqs:
                print >> w_f,'>{0}'.format(pro)
                print >> w_f,'{0}'.format(seq)


if __name__ == "__main__":
    main()
