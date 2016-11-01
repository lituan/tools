"""
use pymol to to pair_wise align structures
put the PDBs in a directory
after launching pymol, go to the directory
then source this script
then run 'align_all . '
"""

import sys
import os
from pymol import cmd


def align_all(directory):

    def get_names(directory):

        def files_in_dir(directory):
            for root, dirs, files in os.walk(directory):
                for f in files:
                    if '.pdb' in f or '.PDB' in f:
                        yield os.path.join(root, f)

        names = []
        for f in files_in_dir(directory):
            f_path,f_name = os.path.split(f)
            f_name,f_extention = os.path.splitext(f_name)
            names.append(f_name)
        return names

    names = get_names(directory)

    for n in names:
        print n

    aligns = []
    for i in range(len(names)):
        for j in range(i):
            pro1 = names[i]
            pro2 = names[j]
            cmd.delete('all')
            cmd.load('%s.pdb'%pro1)
            cmd.load('%s.pdb'%pro2)
            # align = cmd.super(pro1,pro2,cycles=20)
            # align = cmd.cealign(pro1,pro2)
            # only align backbone
            align = cmd.super(pro1+' and bb.',pro2+' and bb.')
            aligns.append((pro1,pro2,align))

    # with open('super_align_result.txt','w') as w_f:
    with open('backbone_super_result.txt','w') as w_f:
        for pro1,pro2,align in aligns:
            print >> w_f,'{0:<20}{1:<20}{2:<20}'.format(pro1,pro2,align)


cmd.extend('align_all',align_all)

