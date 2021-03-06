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


def filter_pdb(pdbids,cutoff=3.0):
    """
    pdbids be a string '1A0R:1,1B9X:1,1B9Y:1,1C15:1,1CWW:1,1CY5:1'
    return format:
    [['structureId','chainId','uniprotAcc','resolution','chainLength','releaseDate'],
     ['"1A0R"', '"B"', '"P62871"', '"2.8"', '"340"', '"1998-12-30"'],
     ['"1B9X"', '"A"', '"P62871"', '"3.0"', '"340"', '"1999-02-23"']]
    """

    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbids,
        'customReportColumns':'structureId,uniprotAcc,entityId,resolution,chainLength,releaseDate',
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
    lines = [line for line in lines[1:] if line[4]]
    pdbids = [line[0] for line in lines if float(line[4]) <= cutoff]

    return pdbids

def pdb_method(pdbids):
    """
    pdbids be a string '1A0R:1,1B9X:1,1B9Y:1,1C15:1,1CWW:1,1CY5:1'
    return format:
    [['structureId','chainId','uniprotAcc','resolution','chainLength','releaseDate'],
     ['"1A0R"', '"B"', '"P62871"', '"2.8"', '"340"', '"1998-12-30"'],
     ['"1B9X"', '"A"', '"P62871"', '"3.0"', '"340"', '"1999-02-23"']]
    """

    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbids,
        'customReportColumns':'structureId,resolution,experimentalTechnique',
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
    # lines = [line for line in lines[1:] if line[1]]
    # pdbids = [line[0] for line in lines if float(line[4]) <= cutoff]

    pdbids = [(line[0],line[2]) for line in lines[1:]]

    return pdbids

def rcsb_phos(p):
    """
    quey RCSB for structures containing phosphorylated residues
    """
    url = 'http://www.rcsb.org/pdb/rest/search'

    phos_type,poly_type = p
    query_head = "<orgPdbQuery><version>head</version><queryType>org.pdb.query.simple.ChemCompIdQuery</queryType>"
    query_phos = "<description>Chemical ID(s):  "+phos_type+" and Polymeric type is "+poly_type+"</description>"+"<chemCompId>"+phos_type+"</chemCompId>"
    query_poly = "<polymericType>"+poly_type+"</polymericType></orgPdbQuery>"

    query = query_head + query_phos + query_poly
    req = urllib2.Request(url,data=query)
    response = urllib2.urlopen(req)
    pdbids = response.read()
    pdbids = pdbids.replace('\n',',')
    return (phos_type,poly_type,pdbids)


def main():
    phos_types = ['PTR','TPO','SEP']
    poly_types = ['Polymeric']
    parameters = [(phos,poly) for phos in phos_types for poly in poly_types]

    p = Pool(6)
    r = p.map(rcsb_phos,parameters)
    p.close


    method_pdbs = {}
    for phosi,polyi,pdbids in r:
        if pdbids:
            pdbids = pdb_method(pdbids)
            for p in pdbids:
                if not p[1] in method_pdbs.keys():
                    method_pdbs[p[1]] = [p[0]]
                else:
                    method_pdbs[p[1]].append(p[0])


    methods = []
    pdbs = []
    for method,pdb in method_pdbs.iteritems():
        print method,len(pdb)
        methods.append(method)
        pdbs.append(pdb)
        with open(method+'.txt','w') as w_f:
            for p in pdb:
                print >> w_f,p


    df = pd.DataFrame({'Methods':methods,'PDB Num':map(len,pdbs)})
    f,ax = plt.subplots()
    sns.set_color_codes('pastel')
    sns.barplot(x='Methods',y='PDB Num',data=df,color='b')
    ax.set(title='Num of PDBs using Different Methods')
    # plt.xticks(rotation=90)
    plt.savefig('Num of PDBs using Different Methods',dpi=300)
    plt.close('all')


if __name__ == "__main__":
    main()
