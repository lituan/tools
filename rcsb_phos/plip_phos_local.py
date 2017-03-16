#! /usr/bin/env python
"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
plipcmd - Main script for PLIP command line execution.
Copyright 2014-2015 Sebastian Salentin

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

modified by Li Tuan
"""

# Compatibility
from __future__ import print_function

# Own modules
try:
    from plip.modules.preparation import *
    from plip.modules.visualize import visualize_in_pymol
    from plip.modules.plipremote import VisualizerData
    from plip.modules.report import StructureReport,__version__
    from plip.modules import config
    from plip.modules.mp import parallel_fn
    from plip.modules.webservices import check_pdb_status, fetch_pdb
except ImportError:
    from modules.preparation import *
    from modules.visualize import visualize_in_pymol
    from modules.plipremote import VisualizerData
    from modules.report import StructureReport, __version__
    from modules import config
    from modules.mp import parallel_fn
    from modules.webservices import check_pdb_status, fetch_pdb

# Python standard library
import sys
import argparse
from argparse import ArgumentParser
import time
import multiprocessing
import json
# External libraries
import lxml.etree as et

descript = "Protein-Ligand Interaction Profiler (PLIP) v%s " \
           "is a command-line based tool to analyze interactions in a protein-ligand complex. " \
           "If you are using PLIP in your work, please cite: " \
           "Salentin,S. et al. PLIP: fully automated protein-ligand interaction profiler. " \
           "Nucl. Acids Res. (1 July 2015) 43 (W1): W443-W447. doi: 10.1093/nar/gkv315" % __version__


def threshold_limiter(aparser, arg):
    arg = float(arg)
    if arg <= 0:
        aparser.error("All thresholds have to be values larger than zero.")
    return arg


def process_pdb(pdbfile, outpath):
    write_message(outpath)
    """Analysis of a single PDB file. Can generate textual reports XML, PyMOL session files and images as output."""
    startmessage = '\nStarting analysis of %s\n' % pdbfile.split('/')[-1]
    write_message(startmessage)
    write_message('='*len(startmessage)+'\n')
    mol = PDBComplex()
    mol.output_path = outpath
    mol.load_pdb(pdbfile)
    # #@todo Offers possibility for filter function from command line (by ligand chain, position, hetid)
    for ligand in mol.ligands:
        mol.characterize_complex(ligand)

    create_folder_if_not_exists(outpath)

    # Generate the report files
    streport = StructureReport(mol)

    config.MAXTHREADS = min(config.MAXTHREADS, len(mol.interaction_sets))


    ######################################
    # PyMOL Visualization (parallelized) #
    ######################################

    if config.PYMOL or config.PICS:
        complexes = [VisualizerData(mol, site) for site in sorted(mol.interaction_sets)
                     if not len(mol.interaction_sets[site].interacting_res) == 0]
        if config.MAXTHREADS > 1:
            write_message('\nGenerating visualizations in parallel on %i cores ...' % config.MAXTHREADS)
            parfn = parallel_fn(visualize_in_pymol)
            parfn(complexes, processes=config.MAXTHREADS)
        else:
            [visualize_in_pymol(plcomplex) for plcomplex in complexes]

    if config.XML:  # Generate report in xml format
        streport.write_xml()

    if config.TXT:  # Generate report in txt (rst) format
        streport.write_txt()

def download_structure(inputpdbid):
    """Given a PDB ID, downloads the corresponding PDB structure.
    Checks for validity of ID and handles error while downloading.
    Returns the path of the downloaded file."""
    try:
        if len(inputpdbid) != 4 or extract_pdbid(inputpdbid.lower()) == 'UnknownProtein':
            sysexit(3, 'Invalid PDB ID (Wrong format)\n')
        pdbfile, pdbid = fetch_pdb(inputpdbid.lower())
        pdbpath = tilde_expansion('%s/%s.pdb' % (config.BASEPATH.rstrip('/'), pdbid))
        create_folder_if_not_exists(config.BASEPATH)
        with open(pdbpath, 'w') as g:
            g.write(pdbfile)
        write_message('file downloaded as %s\n\n' % pdbpath)
        return pdbpath, pdbid

    except ValueError:  # Invalid PDB ID, cannot fetch from RCBS server
        sysexit(3, 'Invalid PDB ID (Entry does not exist)\n')

def remove_duplicates(slist):
    """Checks input lists for duplicates and returns
    a list with unique entries"""
    unique = list(set(slist))
    difference = len(slist) - len(unique)
    if difference == 1:
        write_message("Removed one duplicate entry from input list.\n")
    if difference > 1:
        write_message("Removed %i duplicate entries from input list.\n" % difference)
    return unique


def plip_main(inputstructs, inputpdbids):
    """Main function. Calls functions for processing, report generation and visualization."""
    pdbid, pdbpath = None, None
    # #@todo For multiprocessing, implement better stacktracing for errors
    # Print title and version
    title = "* Protein-Ligand Interaction Profiler v%s *" % __version__
    write_message('\n' + '*' * len(title) + '\n')
    write_message(title)
    write_message('\n' + '*' * len(title) + '\n\n')

    if inputstructs is not None:  # Process PDB file(s)
        num_structures = len(inputstructs)
        inputstructs = remove_duplicates(inputstructs)
        for inputstruct in inputstructs:
            if os.path.getsize(inputstruct) == 0:
                sysexit(2, 'Empty PDB file\n')  # Exit if input file is empty
            if num_structures > 1:
                basename = inputstruct.split('.')[0].split('/')[-1]
                config.OUTPATH = '/'.join([config.BASEPATH, basename])
            process_pdb(inputstruct, config.OUTPATH)
    else:  # Try to fetch the current PDB structure(s) directly from the RCBS server
        num_pdbids = len(inputpdbids)
        inputpdbids =remove_duplicates(inputpdbids)
        for inputpdbid in inputpdbids:
            pdbpath, pdbid = download_structure(inputpdbid)
            if num_pdbids > 1:
                # config.OUTPATH = '/'.join([config.BASEPATH, pdbid[1:3].upper(), pdbid.upper()])
                config.OUTPATH = '/'.join([config.BASEPATH, pdbid.upper()])
            process_pdb(pdbpath, config.OUTPATH)

    if (pdbid is not None or inputstructs is not None) and config.BASEPATH is not None:
        if config.BASEPATH in ['.', './']:
            write_message('\nFinished analysis. Find the result files in the working directory.\n\n')
        else:
            write_message('\nFinished analysis. Find the result files in %s\n\n' % config.BASEPATH)

# if __name__ == '__main__':

    ##############################
    # Parse command line arguments
    ##############################

def myplip(pdbf):
    # pdbids format '1NEX,1ARJ...'
    file_path,file_name = os.path.split(pdbf)
    pathname = file_path+'/'+file_name.split('.')[0].upper()
    config.VERBOSE = True
    config.DEBUG = False
    config.MAXTHREADS = 1
    config.XML = True
    config.TXT = True
    config.PICS = True
    config.PYMOL = True
    config.OUTPATH = pathname
    config.OUTPATH = tilde_expansion("".join([config.OUTPATH, '/'])
                                     if not config.OUTPATH.endswith('/') else config.OUTPATH)
    config.BASEPATH = config.OUTPATH  # Used for batch processing
    # config.BASEPATH = 'basetest'
    # config.BREAKCOMPOSITE = arguments.breakcomposite
    # config.ALTLOC = arguments.altlocation
    # config.PEPTIDES = False
    # config.INTRA = arguments.intra
    # config.NOFIX = arguments.nofix
    config.KEEPMOD = True
    # expanded_path = tilde_expansion(arguments.input) if arguments.input is not None else None
    # plip_main(expanded_path, arguments.pdbid)  # Start main script
    plip_main([pdbf], None)  # Start main script,pdbids is a list


if __name__ == "__main__":

    for root,dirs,files in os.walk(sys.argv[-1]):
        for f in files:
            if f[-4:] == '.pdb':
                f = os.path.join(root,f)
                myplip(f)


