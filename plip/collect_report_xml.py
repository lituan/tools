#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
parse report.xml from plip
result has following forms
pdb_lig_interactions format: 1A07 C_101_PTR ['hydrophobic',['A_206_LYS'...]...]

usage: python collect_report_xml.py plip_result_dir
"""

import os
import sys
from bs4 import BeautifulSoup
from collections import OrderedDict
from collections import Counter
import cPickle as pickle

def parse_report(report_f):
    soup = BeautifulSoup(open(report_f),'lxml')

    pdbid = soup.find('pdbid').get_text()
    lig_interactions = OrderedDict()
    for bsite in soup.find_all('bindingsite'):
        lig = bsite.find('chain').get_text() + '_' + bsite.find('position').get_text() + '_' + bsite.find('hetid').get_text()
        lig_interactions[lig] = []

        interactions = bsite.find('interactions')
        hydrophobic = interactions.find_all('hydrophobic_interaction')
        if hydrophobic:
            hydrophobic_res = []
            for r in hydrophobic:
                resn = r.find('resnr').get_text()
                rest = r.find('restype').get_text()
                reschain = r.find('reschain').get_text()
                res = reschain+'_'+resn+'_'+rest
                hydrophobic_res.append(res)
            if hydrophobic_res:
                lig_interactions[lig].append(('hydrophobic',hydrophobic_res))

        hydrogen_bonds = interactions.find_all('hydrogen_bond')
        if hydrogen_bonds:
            hydrogen_bonds_res = []
            for r in hydrogen_bonds:
                resn = r.find('resnr').get_text()
                rest = r.find('restype').get_text()
                reschain = r.find('reschain').get_text()
                res = reschain+'_'+resn+'_'+rest
                hydrogen_bonds_res.append(res)
            if hydrogen_bonds_res:
                lig_interactions[lig].append(('hydrogen_bonds',hydrogen_bonds_res))

        water_bridges = interactions.find_all('water_bridge')
        if water_bridges:
            water_bridges_res = []
            for r in water_bridges:
                resn = r.find('resnr').get_text()
                rest = r.find('restype').get_text()
                reschain = r.find('reschain').get_text()
                res = reschain+'_'+resn+'_'+rest
                water_bridges_res.append(res)
            if water_bridges_res:
                lig_interactions[lig].append(('water_bridges',water_bridges_res))

        salt_bridges = interactions.find_all('salt_bridge')
        if salt_bridges:
            salt_bridges_res = []
            for r in salt_bridges:
                resn = r.find('resnr').get_text()
                rest = r.find('restype').get_text()
                reschain = r.find('reschain').get_text()
                dist =  r.find('dist').get_text()
                res = reschain+'_'+resn+'_'+rest
                salt_bridges_res.append(res)
            if salt_bridges_res:
                lig_interactions[lig].append(('salt_bridges',salt_bridges_res))

        pi_stacks = interactions.find_all('pi_stack')
        if pi_stacks:
            pi_stacks_res = []
            for r in pi_stacks:
                resn = r.find('resnr').get_text()
                rest = r.find('restype').get_text()
                reschain = r.find('reschain').get_text()
                res = reschain+'_'+resn+'_'+rest
                pi_stacks_res.append(res)
            if pi_stacks_res:
                lig_interactions[lig].append(('pi_stacks',pi_stacks_res))

        pi_cation_interactions = interactions.find_all('pi_cation_interaction')
        if pi_cation_interactions:
            pi_cation_interactions_res = []
            for r in pi_cation_interactions:
                resn = r.find('resnr').get_text()
                rest = r.find('restype').get_text()
                reschain = r.find('reschain').get_text()
                res = reschain+'_'+resn+'_'+rest
                pi_cation_interactions_res.append(res)
            if pi_cation_interactions_res:
                lig_interactions[lig].append(('pi_cation_interactions',pi_cation_interactions_res))

        halogen_bonds = interactions.find_all('halogen_bond')
        if halogen_bonds:
            halogen_bonds_res = []
            for r in halogen_bonds:
                resn = r.find('resnr').get_text()
                rest = r.find('restype').get_text()
                reschain = r.find('reschain').get_text()
                res = reschain+'_'+resn+'_'+rest
                halogen_bonds_res.append(res)
            if halogen_bonds_res:
                lig_interactions[lig].append(('halogen_bonds',halogen_bonds_res))

        metal_complexes = interactions.find_all('metal_complex')
        if metal_complexes:
            metal_complexes_res = []
            for r in metal_complexes:
                resn = r.find('resnr').get_text()
                rest = r.find('restype').get_text()
                reschain = r.find('reschain').get_text()
                res = reschain+'_'+resn+'_'+rest
                metal_complexes_res.append(res)
            if metal_complexes_res:
                lig_interactions[lig].append(('metal_complexes',metal_complexes_res))

    # pdb_lig_interactions format: 1A07 C_101_PTR ['hydrophobic',['A_206_LYS'...]...]
    for lig,interactions in lig_interactions.iteritems():
        if not interactions:
            lig_interactions.pop(lig)
    if lig_interactions:
        pdb_lig_interactions = [(pdbid,lig,interactions) for lig,interactions in lig_interactions.iteritems()]
        return pdb_lig_interactions
    else:
        return ''


def main():
    pdb_lig_interactions = []
    for root,dirs,files in os.walk(sys.argv[-1]):
        for f in files:
            if f == 'report.xml':
                fname = os.path.split(root)[-1]
                f = os.path.join(root,f)
                result = parse_report(f)
                if result:
                    for pdb,lig,interactions in result:
                        if 'TPO' in lig or 'SEP' in lig or 'PTR' in lig:
                            pdb_lig_interactions.append([pdb,lig,interactions])

    if pdb_lig_interactions:
        print 'get some pdb_lig interactions'
        fname = '.'.join(os.path.split(sys.argv[-1])[-1].split('.')[:-1])
        pickle.dump(pdb_lig_interactions,open(fname+'_pdb_lig_interactions.pickle','w'))
        with open(fname+'_pdb_lig.txt','w') as w_f:
            for pdb,lig,interactions in pdb_lig_interactions:
                print >> w_f,pdb,lig
                for i in interactions:
                    print >> w_f,i
    else:
        print 'get none pdb_lig interactions'

if __name__ == "__main__":
    main()
