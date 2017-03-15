#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
check strange ligands
"""

import os
import sys
import urllib
import urllib2
import cPickle as pickle
from multiprocessing import Pool


def main():

    pdb_interfaces = pickle.load(open(sys.argv[-1]))

    strange_pdb_interfaces = [p for p in pdb_interfaces if len(p[7][0][0]) != 1 or len(p[7][1][0]) != 1]


    with open('strange_pdb_interfaces.txt','w') as w_f:
        for i in strange_pdb_interfaces:
            pdbid,interface_id,interface_area,inter_solvent,pvalue,binding_energy,interaction_res,chains = i
            print >> w_f, '{0:<8}{1:<8}{2:<8.2f}{3:<8.2f}{4:<8.2f}{5:<8.2f}\t{6:<}\t{7:<}'.format(pdbid,interface_id,interface_area,inter_solvent,pvalue,binding_energy,interaction_res,chains)
    pickle.dump(pdb_interfaces,open('pdb_interfaces.pickle','w'))


if __name__ == "__main__":
    main()

