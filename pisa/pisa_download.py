#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
get interface xml from ebi
"""

import os
import sys
import urllib
import urllib2
import httplib
from collections import Counter
from multiprocessing import Pool


def filter_pdb(pdbids,cutoff=3.0):
    """
    pdbids be a string '1A0R:1,1B9X:1,1B9Y:1,1C15:1,1CWW:1,1CY5:1'
    return format:
    [['structureId','chainId','uniprotAcc','resolution','chainLength','releaseDate'],
     ['"1A0R"', '"B"', '"P62871"', '"2.8"', '"340"', '"1998-12-30"'],
     ['"1B9X"', '"A"', '"P62871"', '"3.0"', '"340"', '"1999-02-23"']]
    """

    pdbids_str = ','.join(pdbids)
    url = 'http://www.rcsb.org/pdb/rest/customReport.csv?'
    data = {
        'pdbids':pdbids_str,
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
    pdbids = pdbids.split(',')[:-1]
    return (phos_type,poly_type,pdbids)

def check_current_directory(directory):
    current_got_pdbids = []
    for root,dirs,files in os.walk(directory):
        for f in files:
            if f[-4:] == '.xml' and len(f.split('_')) == 2:
                current_got_pdbids.append(f[:4].upper())
    return set(current_got_pdbids)

def get_pisa(pdbid):
    url = 'http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?'+pdbid.lower()
    for i in range(10):
        try:
            pdbinterface = urllib2.urlopen(url).readlines()
            with open(pdbid+'_interface.xml','w') as w_f:
                for line in pdbinterface:
                    print >> w_f,line.rstrip('\r\n')
            print pdbid+' download complete'
        except httplib.IncompleteRead, e:
            print pdbid,e
            continue


def main():
    phos_types = ['PTR','TPO','SEP']
    poly_types = ['Polymeric']
    parameters = [(phos,poly) for phos in phos_types for poly in poly_types]

    p = Pool(6)
    r = p.map(rcsb_phos,parameters)
    p.close

    current_got_pdbids = check_current_directory('.')

    resolution_cutoff = 3.0
    parameters = []
    for phosi,polyi,pdbids in r:
        if pdbids:
            pdbids = [p.split(':')[0] for p in pdbids]
            pdbids = filter_pdb(pdbids,resolution_cutoff)
            if pdbids:
                pdbids = [p for p in pdbids if not p in current_got_pdbids]
                if pdbids:
                    parameters += pdbids

    parameters = set(parameters)
    print parameters

    for i in range(20):
        try:
            p = Pool(6)
            r = p.map(get_pisa,parameters)
            p.close()
            break
        except:
            current_got_pdbids = check_current_directory('.')
            parameters = [p for p in parameters if not p in current_got_pdbids]
            p = Pool(6)
            r = p.map(get_pisa,parameters)
            p.close()



if __name__ == "__main__":
    main()
