#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import os
import urllib
import urllib2
import httplib
import lxml.etree as et

def check_pdb_status(pdbids):
    """Returns the status and up-to-date entry in the PDB for a given PDB ID"""
    pdbids = [p.upper() for p in pdbids]
    pdbids = ','.join(pdbids)
    for i in range(10):
        try:
            url = 'http://www.rcsb.org/pdb/rest/idStatus?structureId=%s' % pdbids
            xmlf = urllib2.urlopen(url)
            xml = et.parse(xmlf)
            xmlf.close()

            status = ''
            replaceby = ''
            results = []
            for df in xml.xpath('//record'):
                status = df.attrib['status']  # Status of an entry can be either 'UNKWOWN', 'OBSOLETE', or 'CURRENT'
                if status == 'CURRENT':
                    results.append([df.attrib['structureId'],df.attrib['status']])
                if status == 'OBSOLETE':
                    results.append([df.attrib['structureId'],df.attrib['status'],df.attrib['replaceBy']])
                if status == 'UNKNOWN':
                    results.append([df.attrib['structureId'],df.attrib['status']])

            current_pdbids = [r[0] for r in results if r[1] == 'CURRENT']
            replace_pdbids = [r[2] for r in results if r[1] == 'OBSOLETE']

            obsolete = [r for r in results if r[1] == 'OBSOLETE']
            unknown = [r for r in results if r[1] == 'UNKNOWN']
            if obsolete:
                with open('obsolete.txt','w') as w_f:
                    for r in obsolete:
                        print >> w_f,'{0:<10}{1:<10}{2:<10}'.format('pdbid','status','replaceby')
                        print >> w_f,'{0:<10}{1:<10}{2:<10}'.format(r[0],r[1],r[2])
            if unknown:
                with open('unknown.txt','w') as w_f:
                    for r in unknown:
                        print >> w_f,'unknown pdb ids'
                        for p in unknown:
                            print >> w_f,p

            return current_pdbids + replace_pdbids

        except Exception,e:
            print e

    print 'something is wrong,please check scripts'
    return 0


def check_current_directory(pdbids,pdbpath):
    if not os.path.exists(pdbpath):
        return pdbids
    current = []
    for root,dirs,files in os.walk(pdbpath):
        for f in files:
            if f[-4:] == '.pdb' and len(f) == 8:
                current.append(f[:-4].upper())
    pdbids = [p for p in pdbids if not p.upper() in current]
    return pdbids


def fetch_pdb(p):
    pdbid,pdbpath = p
    print 'Downloading file from PDB...'
    pdb_format = 'cif'
    for i in range(10):
        try:
            pdburl = 'http://files.rcsb.org/download/'+pdbid.lower()+'.'+pdb_format
            pdbfile = urllib2.urlopen(pdburl).read()
            if 'sorry' in pdbfile:
                print 'No file in PDB format available from wwwPDB for the given PDB ID ' + pdbid
                return pdbid,'sorry'
            else:
                if not os.path.exists(pdbpath):
                    os.makedirs(pdbpath)
                pdbpath = os.path.join(pdbpath,pdbid+'.'+pdb_format)
                with open(pdbpath,'w') as w_f:
                    w_f.write(pdbfile)
                    print pdbid+' is downloaded'
                    return pdbid,'success'
        except urllib2.HTTPError as e:
            print e
            print 'No file in PDB format available from wwwPDB for the given PDB ID ' + pdbid
            pdb_format = 'cif'

    return pdbid,'fail'


def download_pdbs(pdbids,pdbpath):
    pdbids = [p.upper() for p in pdbids]
    input_num = len(pdbids)
    print 'input pdb ids',input_num

    pdbids = check_current_directory(pdbids,pdbpath)
    rest_num = len(pdbids)
    print input_num - rest_num,'already downlaoded'
    if len(pdbids) == 0:
        return

    pdbids = check_pdb_status(pdbids)
    current_num = len(pdbids)
    if len(pdbids) == 0:
        print 0,'to be downlaod'
        return
    print current_num,'is going to be downloaded'

    for i in range(10):
        parameters = [(p,pdbpath) for p in pdbids]
        # p = Pool(8)
        # result = p.map(fetch_pdb,parameters)
        # p.close()
        result = []
        for p in parameters:
            r = fetch_pdb(p)
            result.append(r)

        fail = [r[0] for r in result if r[1] == 'fail']
        if len(fail) == 0:
            break
        else:
            pdbids = fail

