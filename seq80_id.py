"reform sequence so that each line contain 80 residues"

import os
import sys

def read_id(id_f):
    with open(id_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line for line in lines if line]
        lines = [line for line in lines]

        return lines

def readfa(fa_f):
    # readin seqs in fasta format
    # seqs foramt:[(pro,seq),...]
    lines = fa_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    pro_line_num = [i for i, line in enumerate(lines) if '>' in line] + [len(lines)]
    seqs = [lines[n:pro_line_num[i + 1]]
            for i, n in enumerate(pro_line_num[:-1])]
    seqs = [(seq[0][1:], ''.join(seq[1:])) for seq in seqs]
    return seqs

def writefa(seqs,ids):

    with open('seq80.fa','w') as w_f:

        # seqs = [(pro.split()[0].split('|')[1],seq) for pro,seq in seqs]
        # seqs = sorted(seqs)
        for i in ids:
            for pro,seq in seqs:
                if pro == i:
                    print >> w_f,'>{}'.format(pro)
                    for i in [seq[i:i+80] for i in range(0,len(seq),80)]:
                        print >> w_f,i
ids = read_id(sys.argv[-2])
with open(sys.argv[-1]) as fa_f:
    seqs = readfa(fa_f)
    writefa(seqs,ids)
