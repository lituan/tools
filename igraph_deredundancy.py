#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import lt

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

    def pim(seq1,seq2):
        identity = len([i for i,s in enumerate(seq1) if s == seq2[i]])
        return identity*1.0/len(seq1)

    scores = []
    seqlen = len(seqs)
    for i in range(seqlen):
        score_i = []
        for j in range(seqlen):
            if j < i:
                score_i.append(scores[j][i])
            elif j > i:
                score_i.append(pim(seqs[i][1],seqs[j][1]))
            else:
                score_i.append(1.0)
        scores.append(score_i)
    return scores

@lt.run_time
def igraph_mis(seqs,scores,cutoff,good_list=[],bad_list=[],no_bad=False):
    """
    use igraph to find maximal independent sets
    sort mis according to quality
    no_bad=False, choose the one with least bad entries
        quality score: length + good_number*0.1-bad_number*0.1
    no_bad=True, no bad entry is allowed
        quality score: length + good_number*0.1-bad_number
    """
    import igraph
    # score greater than cutoff indicates a link
    adj_m = [map(lambda x: 1 if x > cutoff else 0,row) for row in scores]
    # no circle
    for i in range(len(adj_m)):
        adj_m[i][i] = 0
    graph = igraph.Graph.Adjacency(adj_m)
    mis = graph.independent_vertex_sets()

    total = len(labels)
    good_list_index = set([labels.index(g) for g in good_list])
    bad_list_index = set([labels.index(b) for b in bad_list])

    if no_bad == False:
        mis_quality = [len(m)+len(set(m).intersection(good_list_index))*1.0/total-len(set(m).intersection(bad_list_index))*1.0/total for m in mis]
        mis = sorted(mis,key=lambda x: mis_quality[mis.index(x)],reverse=True)
    elif no_bad == True:
        mis_quality = [len(m)+len(set(m).intersection(good_list_index))*1.0/total-len(set(m).intersection(bad_list_index)) for m in mis]
        mis = sorted(mis,key=lambda x: mis_quality[mis.index(x)],reverse=True)

    mmis = mis[0]

    nr_seqs = [seqs[i] for i in mmis ]
    nr_scores = [[score_i[i] for i in mmis] for score_i in scores]
    nr_scores = [nr_scores[i] for i in mmis]
    return nr_seqs,nr_scores

def plot_heatmap(seqs, scores,file_name):

    if not len(scores) > 2:
        print 'number of seqs is too small'
        return

    column_labels = [s[0] for s in seqs]
    row_labels = column_labels
    scores = [map(lambda x: float(x), row) for row in scores]
    scores = np.array(scores)
    distances = [map(lambda x: 1-x,row) for row in scores]
    linkage = sch.linkage(spd.squareform(distances),method='average')
    df = pd.DataFrame(scores,columns=column_labels, index=row_labels)

    if len(df.columns) > 30:
        # sns_plot = sns.clustermap(df,row_linkage=linkage,col_linkage=linkage,annot=True,fmt='3.2f',xticklabels='',yticklabels='')
        sns_plot = sns.clustermap(df,row_linkage=linkage,col_linkage=linkage,xticklabels='',yticklabels='')
    else:
        sns_plot = sns.clustermap(df,row_linkage=linkage,col_linkage=linkage,annot=True,fmt='3.2f')
        # sns_plot = sns.clustermap(df,row_linkage=linkage,col_linkage=linkage)
        plt.setp(sns_plot.ax_heatmap.yaxis.get_majorticklabels(), rotation=20)
        plt.setp(sns_plot.ax_heatmap.xaxis.get_majorticklabels(), rotation=70)
    # plt.yticks(rotation=90)
    sns_plot.savefig(file_name+'.png',dpi=300)
    plt.close('all')

def main():
    cutoff = 0.9
    seqs = read_msa(sys.argv[-1])

    scores = get_pim(seqs)
    nr_seqs,nr_scores = igraph_mis(seqs,scores,cutoff)

    plot_heatmap(nr_seqs,nr_scores,'nr_seqs_'+'_'+str(cutoff))

    # write_msa(nr_seqs,'nr_seqs_'+str(cutoff))

main()






