"""
transform percent identify matrix from clustal into other forms

"""

import os
import sys

def read_pim(pim_f):
    lines = pim_f.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line for line in lines if len(line) > 0]

    scores = []
    ids = []

    ids = [line.split()[1] for line in lines if ':' in line]
    scores = [[float(i) for i in line.split()[2:]] for line in lines if ':' in line]

    return ids,scores

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

def write_resutls(ids, scores,file_path,file_name):
    result = [[i] + score for i, score in zip(ids, scores)]
    header = [['ID'] + ids ]
    result = header + result

    filename = os.path.join(file_path,file_name+'_scores_tab.txt')
    with open(filename,'w') as w_f:
        for r in result:
            print >> w_f, '\t'.join([str(ri)for ri in r])

    result = align_lis_lis(result)
    filename = os.path.join(file_path,file_name+'_scores_align.txt')
    with open(filename,'w') as w_f:
        for r in result:
            print >> w_f, '\t'.join([str(ri)for ri in r])

    pair_scores = []
    for i in range(len(ids)):
        for j in range(i):
            pair_scores.append([scores[i][j],(ids[i],ids[j])])

    pair_scores = sorted(pair_scores)
    filename = os.path.join(file_path,file_name+'_pair_scores.txt')
    with open(filename,'w') as w_f:
        for score,pair in pair_scores:
            print >> w_f,score,'\t','{0:<25}{1:<}'.format(pair[0],pair[1])

def main():

    f = sys.argv[-1]
    f_path,f_full_name = os.path.split(f)
    f_name,f_exten = os.path.splitext(f_full_name)

    with open(f) as pim_f:
        ids,scores = read_pim(pim_f)
        write_resutls(ids,scores,f_path,f_name)


if __name__ == "__main__":
    main()
