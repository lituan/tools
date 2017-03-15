#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run pisa to parse interfaces in PDBs
first install CCP4 and change the following variables
"""
import os
import sys
import subprocess
import cPickle as pickle
from multiprocessing import Pool

os.environ['CCP4_SCR'] = 'C:\\ccp4temp'
os.environ['CCP4I_TCLTK'] = 'C:\\CCP4-7\\TclTk84\\bin'
os.environ['CBIN'] = 'C:\\CCP4-7\\7.0\\bin'
os.environ['CLIB'] = 'C:\\CCP4-7\\lib'
os.environ['CLIBD'] = 'C:\\CCP4-7\\lib\\data'
os.environ['CEXAM'] = 'C:\\CCP4-7\\examples'
os.environ['CHTML'] = 'C:\\CCP4-7\\html'
os.environ['CINCL'] = 'C:\\CCP4-7\\include'
os.environ['CCP4I_TOP'] = 'C:\\CCP4-7\\share\\ccp4i'
os.environ['CLIBD_MON'] = 'C:\\CCP4-7\\lib\\data\\monomers\\'
os.environ['MMCIFDIC'] = 'C:\\CCP4-7\\lib\\ccp4\\cif_mmdic.lib'
os.environ['CRANK'] = 'C:\\CCP4-7\\share\\ccp4i\\crank'
os.environ['CCP4_OPEN'] = 'unknown'
os.environ['GFORTRAN_UNBUFFERED_PRECONNECTED'] = 'Y'

os.environ['PATH'] = 'C:\\CCP4-7\\7.0\\bin'
os.environ['PISA_CONF_FILE'] = 'C:\\CCP4-7\\7.0\\share\\pisa\\pisa.cfg'

def pisa(f):
    if not os.path.exists('detail'):
        os.makedirs('detail')
    if not os.path.exists('interface_xml'):
        os.makedirs('interface_xml')
    pdbid = f[-8:-4].lower()
    # subprocess.call(['pisa',pdbid,'-analyse',f])
    # xml_fname = os.path.join('interface_xml',pdbid+'_inteface.xml')
    # subprocess.call(['pisa',pdbid,'-xml','interfaces','>',xml_fname],shell=True)
        # output = subprocess.check_output(['pisa',pdbid,'-detail','interfaces',str(interface_num)],shell=True)

    for interface_num in range(100,200):
        try:
            output = subprocess.check_output(['pisa',pdbid,'-detail','interfaces',str(interface_num)],shell=True)

            detail_fname = os.path.join('detail',pdbid+'_'+str(interface_num)+'_detail.txt')
            subprocess.call(['pisa',pdbid,'-detail','interfaces',str(interface_num),'>',detail_fname],shell=True)
        except:
            continue


def main():
    parameters = []
    for root,dirs,files in os.walk(sys.argv[-1]):
        for f in files:
            if f[-4:] == '.pdb' and len(f) == 8:
                f = os.path.join(root,f)
                parameters.append(f)

    p = Pool(8)
    p.map(pisa,parameters)
    p.close()

if __name__ == "__main__":
    main()


