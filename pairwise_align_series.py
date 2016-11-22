import os
import sys
from wdsp import Wdsp

"""
input seqs in fasta format, compare the rest sequences with the firs one,
rank the sequences according to the similarity,
output sequence similarity series
usuage example:
python pairwise_align example.fasta
"""

def align_lis_lis(lis_lis):
    """align and  nested list to print a table"""
    lis_lis = [[str(l) for l in lis]
               for lis in lis_lis]  # trans every element to str
    #make all inner lists of the same length
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [lis + (inner_lis_max_len - len(lis)) * [''] for lis in lis_lis]
    #trans list, so that the elements of the same column are in one list
    lis_lis = [[lis[i] for lis in lis_lis] for i in range(inner_lis_max_len)]
    #make element in the same list have the same length
    aligned = []
    for lis in lis_lis:
        width = max([len(l) for l in lis])
        lis = [l + (width - len(l)) * ' ' for l in lis]
        aligned.append(lis)
    #trans list_list to the original list_list
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [[lis[i] for lis in aligned] for i in range(inner_lis_max_len)]
    return lis_lis

def align(seq1, seq2):
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62
    gap_open = -10  # usual value
    gap_extend = -0.5  # usual value

    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)

    seq1 = alns[0][0]
    seq2 = alns[0][1]
    identity = [1 for i, s in enumerate(seq1) if s == seq2[i]]
    identity = 1.0 * len(identity) / len(seq1)

    return '{0:<12.10f}'.format(identity)


def readfa(fa_f):
    # readin seqs in fasta format
    # seqs foramt:[(pro,seq),...]
    lines = fa_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    pro_line_num = [i for i, line in enumerate(
        lines) if '>' in line] + [len(lines)]
    seqs = [lines[n:pro_line_num[i + 1]]
            for i, n in enumerate(pro_line_num[:-1])]
    seqs = [(seq[0][1:], ''.join(seq[1:])) for seq in seqs]
    return seqs

def read_wdsp(wdsp_f):
    w = Wdsp(wdsp_f)
    seqs = w.seqs
    seqs = [(k,v) for k,v in seqs.iteritems()]
    return seqs


def align_seqs(seqs):
    scores = [(align(seqs[0][1],v),k,v) for (k,v) in seqs]
    scores = sorted(scores,reverse=True)
    return scores


def write_resutls(seqs, scores,file_path,file_name):

    filename = os.path.join(file_path,file_name+'scores_series.txt')
    with open(filename,'w') as w_f:
        for s in scores:
            print >> w_f, '\t'.join([str(ri)for ri in [s[1],s[0],s[2]]])

def main():
    with open(sys.argv[-1])  as fa_f:
        seqs = readfa(fa_f)
    # with open(sys.argv[-1]) as wdsp_f:
        # seqs = read_wdsp(wdsp_f)

    scores = align_seqs(seqs)

    file_path,file_name = os.path.split(sys.argv[-1])
    file_name,file_extention = os.path.splitext(file_name)
    write_resutls(seqs, scores,file_path,file_name)
    # plot_heatmap(seqs,scores)

if __name__ == "__main__":
    main()
