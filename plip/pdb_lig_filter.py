#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
find interactions involving only one interacting chain
"""

import os
import sys
import cPickle as pickle
from collections import Counter


def filter_pdb_lig(pdb_lig_interactions):

    pdb_lig_interaction_res= []
    for pdb,lig,interactions in pdb_lig_interactions:
        inter_res = []
        for interaction_type,interaction_res in interactions:
            inter_res += interaction_res
        if inter_res:
            pdb_lig_interaction_res.append((pdb,lig,set(inter_res)))

    pdb_lig_one_chain = []
    for pdb,lig,res in pdb_lig_interaction_res:
        chains = [r.split('_')[0] for r in res]
        chain_num = len(set(chains))
        if chain_num == 1:
            pdb_lig_one_chain.append((pdb,lig))

    pdb_lig_interactions_one_chain = []
    for pdb,lig,interactions in pdb_lig_interactions:
        if (pdb,lig) in pdb_lig_one_chain:
            pdb_lig_interactions_one_chain.append((pdb,lig,interactions))

    print 'one_chain pdb_lig', len(pdb_lig_interactions_one_chain)
    return pdb_lig_interactions_one_chain


def select_interaction(pdb_lig_interactions):
    pdb_lig_interactions_select = []
    for pdb,lig,interactions in pdb_lig_interactions:
        interactions_select = []
        for interaction in interactions:
            if interaction[0] in ['hydrogen_bonds','salt_bridges','water_bridges']:
                interactions_select.append(interaction)
        if interactions_select:
            if pdb == '1OL5':
                print pdb,lig,interactions_select
            pdb_lig_interactions_select.append((pdb,lig,interactions_select))
        if pdb == '1OL5':
            print 'ddd'
            # print pdb,lig,interactions_select
    return pdb_lig_interactions_select

def main():

    pdb_lig_interactions = pickle.load(open(sys.argv[-1],'r'))

    pdb_lig_interactions_select = select_interaction(pdb_lig_interactions)

    pdb_lig_interactions_filter = filter_pdb_lig(pdb_lig_interactions_select)

    fname = '.'.join(os.path.split(sys.argv[-1])[-1].split('.')[:-1])
    pickle.dump(pdb_lig_interactions_filter,open(fname+'_filter'+'.pickle','w'))


if __name__ == "__main__":
    main()
