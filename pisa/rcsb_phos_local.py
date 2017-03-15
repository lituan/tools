#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
use pfam accs and chemical ID to search structures Phos-binding
"""
import os
import sys
import urllib
import urllib2
from multiprocessing import Pool
import lxml.etree as et


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
        'pdbids':','.join(pdbids),
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



def write_lis_lis(lis_lis,filename,cols=[]):
    """align nested list to print a table"""

    lis_lis = [lis if lis else ['    '] for lis in lis_lis]
    lis_lis = [[str(l) for l in lis]
               for lis in lis_lis]  # trans every element to str
    #make all inner lists of the same length
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [lis + (inner_lis_max_len - len(lis)) * [''] for lis in lis_lis]
    #make element in the same list have the same length
    aligned = []
    for lis in lis_lis:
        width = max([len(l) for l in lis])
        lis = [l + (width - len(l)) * ' ' for l in lis]
        aligned.append(lis)
    new_lis_lis = [';'.join([aligned[i][j] for i in range(len(aligned))]) for j in range(len(aligned[0]))]
    with open(filename+'.txt','w') as w_f:
        if cols:
            print >> w_f,'\t;'.join(cols)
        for l in new_lis_lis:
            print >> w_f,l

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
    pdbids = pdbids.split(',')
    pdbids = [p for p in pdbids if p]
    return (phos_type,poly_type,pdbids)

def check_pdb_status(pdbid):
    """Returns the status and up-to-date entry in the PDB for a given PDB ID"""
    url = 'http://www.rcsb.org/pdb/rest/idStatus?structureId=%s' % pdbid
    xmlf = urllib2.urlopen(url)
    xml = et.parse(xmlf)
    xmlf.close()
    status = None
    current_pdbid = pdbid
    for df in xml.xpath('//record'):
        status = df.attrib['status']  # Status of an entry can be either 'UNKWOWN', 'OBSOLETE', or 'CURRENT'
        if status == 'OBSOLETE':
            current_pdbid = df.attrib['replacedBy']  # Contains the up-to-date PDB ID for obsolete entries
    return [status, current_pdbid.lower()]

def fetch_pdb(p):
    pdbid,pdbpath = p
    pdbid = pdbid.lower()
    print 'Checking status of PDB ID ' + pdbid
    state,current_entry = check_pdb_status(pdbid)
    if state == 'OBSOLETE':
        print 'entry is obsolete, getting ' + current_entry + ' instead'
    elif state == 'CURRENT':
        print 'entry is up to date'
    elif state == 'UNKNOWN':
        print 'Invalid PDB ID ' + pdbid
    print 'Downloading file from PDB...'
    pdburl = 'http://www.rcsb.org/pdb/files/'+current_entry+'.pdb'
    try:
        pdbfile = urllib2.urlopen(pdburl).read()
        if 'sorry' in pdbfile:
            print 'No file in PDB format available from wwwPDB for the given PDB ID ' + current_entry
        else:
            if not os.path.exists(pdbpath):
                os.makedirs(pdbpath)
            pdbpath = os.path.join(pdbpath,pdbid+'.pdb')
            with open(pdbpath,'w') as w_f:
                w_f.write(pdbfile)
                print pdbid+' is downloaded'
                return pdbpath
    except urllib2.HTTPError:
        print 'No file in PDB format available from wwwPDB for the given PDB ID ' + current_entry

def check_current_directory(directory):
    current_got_pdbids = []
    for root,dirs,files in os.walk(directory):
        for f in files:
            if f[-4:] == '.pdb' and len(f) == 8:
                current_got_pdbids.append(f[:4].upper())
    return set(current_got_pdbids)

def main():
    phos_types = ['PTR','TPO','SEP']
    poly_types = ['Polymeric']
    parameters = [(phos,poly) for phos in phos_types for poly in poly_types]

    p = Pool(6)
    r = p.map(rcsb_phos,parameters)
    p.close


    for i in range(100):
        current_got_pdbids = check_current_directory('.')

        resolution_cutoff = 3.5
        parameters = []
        for phosi,polyi,pdbids in r:
            if pdbids:
                pdbids = filter_pdb(pdbids,resolution_cutoff)
                if pdbids:
                    for p in pdbids:
                        if not p in current_got_pdbids:
                            if not (p,'pdb') in parameters:
                                parameters.append((p,'pdb'))
        if not parameters:
            break

        p = Pool(8)
        r = p.map(fetch_pdb,parameters)
        p.close()


if __name__ == "__main__":
    main()
