#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
get branch length from a newick tree
"""
from ete3 import Tree

t = Tree('organisms_in_timetree_in_wdsp.nwk',format=1)
leafs = t.get_leaf_names()
hd = []
spe = 'Homo_sapiens'
for l in leafs:
    if not l == spe:
        hd.append((spe,l,t.get_distance(spe,l)))

hd = sorted(hd,key=lambda x: x[2])
with open('distance.txt','w') as w_f:
    for spe,l,d in hd:
        print >> w_f, '{0:<30}{1:<30}{2:<}'.format(spe,l,d)
