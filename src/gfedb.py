#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TODO: **** ADD HAS_FEATURE
#       between SEQUENCE and features
# from seqann.models.annotation import Annotation

from seqann.gfe import GFE
import logging
import argparse
import os
import sys
import urllib.request
from gfedb_utils import build_hla_graph
from constants import *
import pandas as pd

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)


def main():
    """This is run if file is directly executed, but not if imported as
    module. Having this in a separate function  allows importing the file
    into interactive python, and still able to execute the
    function for testing"""
    parser = argparse.ArgumentParser()

    parser.add_argument("-k", "--kir",
                        required=False,
                        help="Bool for KIR",
                        action='store_true')

    parser.add_argument("-a", "--align",
                        required=False,
                        help="Bool for loading alignments",
                        action="store_true")

    parser.add_argument("-d", "--debug",
                        required=False,
                        help="Bool for debugging",
                        action="store_true")

    parser.add_argument("-o", "--out_dir",
                        required=True,
                        help="Output directory",
                        type=str,
                        action="store")

    parser.add_argument("-n", "--number",
                        required=False,
                        help="Number of IMGT/DB releases",
                        default=1,
                        type=int,
                        action="store")

    parser.add_argument("-c", "--count",
                        required=False,
                        help="Number of alleles",
                        default=1,
                        type=int,
                        action="store")

    parser.add_argument("-r", "--release",
                        required=False,
                        help="IMGT/DB release",
                        type=str,
                        action="store")

    parser.add_argument("-v", "--verbose",
                        help="Option for running in verbose",
                        action="store_true")

    parser.add_argument("-l", "--limit",
                        required=False,
                        help="Limit number of records in output",
                        default=None,
                        nargs='?',
                        type=int,
                        action="store")

    args = parser.parse_args()

    logging.info(f'args: {vars(args)}')

    #logging.info(f'args:\n{vars(args)}')

    last_release = pd.read_html(imgt_hla)[0]['Release'][0].replace(".", "")
    #logging.info(f'last_release: {last_release}')

    release_n = args.number
    dbversion = args.release if args.release else last_release
    verbosity = 1
    num_alleles = args.count

    debug = True if args.debug else False
    kir = True if args.kir else False
    align = True if args.align else False
    verbose = True if args.verbose else False
    load_loci = hla_loci + kir_loci if kir else hla_loci

    if debug:
        logging.info("Running in debug mode")
        load_loci = ["HLA-A"]
        kir = False
        verbose = True
        verbosity = 2
        release_n = 1

    # # Get last five IMGT/HLA releases
    # if releases:
    #     dbversions = [db for db in releases.split(",")]
    # else:
    #     dbversions = pd.read_html(imgt_hla)[0]['Release'][0:release_n].tolist()

    gfe_maker = GFE(verbose=verbose, 
        verbosity=verbosity,
        load_features=False, 
        store_features=True,
        loci=load_loci)

    # TO DO: for this for loop into the build.sh script
    #for dbversion in dbversions:
    
    #logging.info(f'\n\nBuilding graph for IMGTHLA version {dbversion[0]}.{dbversion[1:3]}.{dbversion[3]}...')
    #logging.info(f'Limit: {args.limit}')
    build_hla_graph(
        dbversion=dbversion, 
        alignments=align, 
        verbose=verbose,
        limit=args.limit,
        num_alleles=num_alleles,
        gfe_maker=gfe_maker)

    logging.info(f'Finished build for version {dbversion[0]}.{dbversion[1:3]}.{dbversion[3]}')

    return

if __name__ == '__main__':
    main()
