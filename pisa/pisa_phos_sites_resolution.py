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
    pdb_phos_sites = pickle.load(open(sys.argv[-1]))

    pdbids = [p[0] for p in pdb_phos_sites]
    pdb_resolution = get_pdb_resolution(pdbids)

    resolutions = [p[1] for p in pdb_resolution]

    f,ax = plt.subplots(figsize=(6,4))
    sns.set_color_codes('pastel')
    b = sns.distplot(resolutions)
    b.figure.subplots_adjust(top=0.9,bottom=0.13,left=0.13,right=0.9)
    ax.set(xlabel='Resolution',ylabel='PDB Num',title='Distribution of Resolutions of final-select-PDBs')
    # plt.xticks(rotation=90)
    plt.savefig(fname+'_resolution_dist.png',dpi=300)
    plt.close('all')




if __name__ == "__main__":
    main()
