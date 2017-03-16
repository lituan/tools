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

def rcsb_pfam(p):
    """
    use beta_sheet, chain length to make sure we get WD40s other than others
    accs be a list, ['Q969H0,P07834']
    return string format: '1A0R:1,1B9X:1,1B9Y:1,1C15:1'
    """
    url = 'http://www.rcsb.org/pdb/rest/search'
    pfam,pfam_dest,pfam_id,phos_type,poly_type = p

    query1 = """
    <orgPdbCompositeQuery version="1.0">
    <queryRefinement>
    <queryRefinementLevel>0</queryRefinementLevel>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.PfamIdQuery</queryType>
    """
    query_pfam = "<description>Pfam Accession Number "+pfam_dest+"</description><pfamID>"+pfam_id+"</pfamID>"

    query2 = """
    </orgPdbQuery>
    </queryRefinement>
    <queryRefinement>
    <queryRefinementLevel>1</queryRefinementLevel>
    <conjunctionType>and</conjunctionType>
    <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ChemCompIdQuery</queryType>
    """
    query_phos = "<description>Chemical ID(s):  "+phos_type+" and Polymeric type is "+poly_type+"</description>"
    query_poly = "<chemCompId>"+phos_type+"</chemCompId>" + "<polymericType>"+poly_type+"</polymericType>"

    query3 = """
    </orgPdbQuery>
    </queryRefinement>
    </orgPdbCompositeQuery>
    """

    query = query1 + query_pfam + query2 + query_phos + query_poly + query3
    req = urllib2.Request(url,data=query)
    response = urllib2.urlopen(req)
    pdbids = response.read()
    pdbids = pdbids.replace('\n',',')
    return (pfam,phos_type,poly_type,pdbids)

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

    pfams = [['SH2','PF00017 SH2 domain','PF00017'],
             ['14333','PF00244 14-3-3 protein','PF00244'],
             ['BRCT','PF00533 BRCA1 C Terminus (BRCT) domain','PF00533'],
             ['FHA','PF00498 FHA domain','PF00498'],
             ['MH2','PF03166 MH2 domain','PF03166'],
             ['PBD','PF00786 P21-Rho-binding domain','PF00786'],
             ['PTB','PF08416 Phosphotyrosine-binding domain','PF08416'],
             ['WW','PF00397 WW domain','PF00397'],
             ['C2','PF00168 C2 domain','PF00168'],
             ['WD40','PF00400 WD domain, G-beta repeat','PF00400']
             ]

    phos_types = ['PTR','TPO','SEP']
    poly_types = ['Polymeric','Any']
    poly_types = ['Polymeric']

    parameters = [(pfam[0],pfam[1],pfam[2],phos_type,poly_type) for pfam in pfams for phos_type in phos_types for poly_type in poly_types]
    p = Pool(4)
    r = p.map(rcsb_pfam,parameters)
    p.close()

    current_got_pdbids = check_current_directory('.')

    resolution_cutoff = 3.0
    parameters = []
    for pfami,phosi,polyi,pdbids in r:
        if pdbids:
            pdbids = pdbids.split(',')[:-1]
            pdbids = [p.split(':')[0] for p in pdbids]
            pdbids = filter_pdb(','.join(pdbids),resolution_cutoff)
            if pdbids:
                for p in pdbids:
                    if not p in current_got_pdbids:
                        parameters.append((p,pfami+'_'+phosi+'_'+polyi))

    # for p in parameters:
        # fetch_pdb(p)
    p = Pool(4)
    r = p.map(fetch_pdb,parameters)
    p.close()


if __name__ == "__main__":
    main()
