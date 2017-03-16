#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
parse pisa assemble xml file
get stable interfaces
"""

import os
import sys
from bs4 import BeautifulSoup
from collections import OrderedDict
from collections import Counter
import cPickle as pickle
from multiprocessing import Pool

def parse_assemble(pisa_f):
    soup = BeautifulSoup(open(pisa_f),'lxml')

    pdb_assemblies= []
    pdb_entries = soup.find_all('pisa_results')
    if pdb_entries:
        for pdb_entry in pdb_entries:
            pdbid = pdb_entry.find('name').get_text().upper()
            assemblies = pdb_entry.find_all('asm_set')
            asms = []
            for asm in assemblies:
                score = asm.find('score').get_text()
                dg_diss = float(asm.find('diss_energy').get_text())
                dg0 = float(asm.find('diss_energy_0').get_text())
                composition = asm.find('composition').get_text()
                interfaces = []
                interface = asm.find('interfaces')
                infs = interface.find_all('interface')
                for inf in infs:
                    inf_id = inf.find('id').get_text()
                    interfaces.append(inf_id)
                if interfaces:
                    asms.append((score,dg_diss,dg0,composition,interfaces))
            if asms:
                pdb_assemblies.append((pdbid,asms))

    return pdb_assemblies


def main():
    parameters = []
    for root,dirs,files in os.walk(sys.argv[-1]):
        for f in files:
            if f[-12:] == 'assemble.xml':
                f = os.path.join(root,f)
                parameters.append(f)

    p = Pool(6)
    result = p.map(parse_assemble,parameters)
    p.close()

    pdb_assemblies = []
    for r in result:
        if r:
            pdb_assemblies += r

    # pdb_assemblies = []
    # for p in parameters:
        # print p
        # pdb_assemblies += parse_assemble(p)


    stable_pdb_interfaces = {}
    for pdbid,assemblies in pdb_assemblies:
        stable_interfaces = []
        for asm in assemblies:
            if asm[0] == 'This assembly appears to be stable in solution.':
                stable_interfaces += asm[4]
        if stable_interfaces:
            stable_pdb_interfaces[pdbid] = set(stable_interfaces)

    pickle.dump(stable_pdb_interfaces,open('stable_pdb_interfaces.pickle','w'))

    with open('stable_pdb_interfaces.txt','w') as w_f:
        for k,v in stable_pdb_interfaces.iteritems():
            print >> w_f, '{0:<8}{1:<}'.format(k,v)

if __name__ == "__main__":
    main()




