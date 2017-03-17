#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
use pfam accs and chemical ID to search structures Phos-binding
"""
import os
import sys
import urllib
import urllib2
import seaborn as sns
import pandas as pd
import numpy as np
from multiprocessing import Pool
import lxml.etree as et
import matplotlib.pyplot as plt
import cPickle as pickle

def get_pdb_resolution(pdbids):

    pdbids = ','.join(pdbids)
    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbids,
        'customReportColumns':'structureId,resolution',
        'service':'wsfile',
        'format':'csv',
    }
    data = urllib.urlencode(data)
    req = urllib2.Request(url,data)
    response = urllib2.urlopen(req)
    lines = response.readlines()
    lines = [line.rstrip('\r\n') for line in lines]
    lines = [line for line in lines if line]
    lines = [line.split(',') for line in lines]
    lines = [[w.strip('"') for w in line] for line in lines]

    # filter pdbs with lower resolution
    lines = [line for line in lines[1:] if line[1]]
    pdb_resolution = [(line[0],float(line[1])) for line in lines if line[1]]

    print 'pdbs with resolution',len(pdb_resolution)
    return pdb_resolution



def main():


    fname = '.'.join(os.path.split(sys.argv[-1])[1].split('.')[:-1])
    pdb_phos_sites1 = pickle.load(open(sys.argv[-3]))
    pdb_phos_sites2 = pickle.load(open(sys.argv[-2]))
    pdb_phos_sites3 = pickle.load(open(sys.argv[-1]))

    pdbids1 = [p[0] for p in pdb_phos_sites1]
    pdb_resolution1 = get_pdb_resolution(pdbids1)
    resolutions1 = [p[1] for p in pdb_resolution1]

    pdbids2 = [p[0] for p in pdb_phos_sites2]
    pdb_resolution2 = get_pdb_resolution(pdbids2)
    resolutions2 = [p[1] for p in pdb_resolution2]

    pdbids3 = [p[0] for p in pdb_phos_sites3]
    pdb_resolution3 = get_pdb_resolution(pdbids3)
    resolutions3 = [p[1] for p in pdb_resolution3]

    f,ax = plt.subplots(figsize=(6,4))
    sns.set_color_codes('pastel')
    b = sns.distplot(resolutions1,color='r',hist=False)
    b = sns.distplot(resolutions2,color='g',hist=False)
    b = sns.distplot(resolutions3,color='b',hist=False)
    b.figure.subplots_adjust(top=0.9,bottom=0.13,left=0.13,right=0.9)
    ax.set(xlabel='Resolution',ylabel='PDB Num',title='Distribution of Resolutions of final-select-PDBs')
    # plt.xticks(rotation=90)
    plt.savefig(fname+'_resolution_dist_3.png',dpi=300)
    plt.close('all')




if __name__ == "__main__":
    main()
