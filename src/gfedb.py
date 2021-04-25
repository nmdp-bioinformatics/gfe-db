#!/usr/bin/env python
# -*- coding: utf-8 -*-

from seqann.gfe import GFE
import logging
import argparse
import os
import sys
import time
import pandas as pd
from gfedb_utils import build_hla_graph
from constants import *

start = time.time()

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

    parser.add_argument("-o", "--out_dir",
                        required=True,
                        help="Output directory",
                        type=str,
                        action="store")

    # TO DO: add option to specify last n releases
    # parser.add_argument("-n", "--number",
    #                     required=False,
    #                     help="Number of IMGT/DB releases",
    #                     default=1,
    #                     type=int,
    #                     action="store")

    parser.add_argument("-c", "--count",
                        required=False,
                        help="Number of alleles",
                        type=int,
                        action="store")

    parser.add_argument("-r", "--release",
                        required=False,
                        help="IMGT/DB release",
                        type=str,
                        action="store")

    parser.add_argument("-p", "--profile",
                        required=False,
                        help="Enable memory profiling",
                        action='store_true')

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

    logging.debug(f'args: {vars(args)}')

    last_release = pd.read_html(imgt_hla)[0]['Release'][0].replace(".", "")
    #logging.info(f'last_release: {last_release}')

    dbversion = args.release if args.release else last_release
    verbosity = 1
    num_alleles = args.count
    limit = args.limit #min(args.count, args.limit)
    mem_profile = args.profile
    kir = True if args.kir else False
    align = True if args.align else False
    verbose = True if args.verbose else False
    

    # # TO DO: Get last n IMGT/HLA releases
    # if releases:
    #     dbversions = [db for db in releases.split(",")]
    # else:
    #     dbversions = pd.read_html(imgt_hla)[0]['Release'][0:release_n].tolist()

    gfe_maker = GFE(verbose=verbose, 
        verbosity=verbosity,
        load_features=False, 
        store_features=True,
        loci=load_loci)
    
    logging.info(f'***************** Building graph for IMGTHLA version {dbversion[0]}.{dbversion[1:3]}.{dbversion[3]} *****************')
    # logging.info(f'Total alleles: {num_alleles}')

    build_hla_graph(
        dbversion=dbversion, 
        alignments=align, 
        verbose=verbose,
        limit=limit,
        #num_alleles=num_alleles,
        gfe_maker=gfe_maker,
        mem_profile=mem_profile)

    logging.info(f'Finished build for version {dbversion[0]}.{dbversion[1:3]}.{dbversion[3]}')
    end = time.time()
    logging.info(f'************************ Build finished in {round(end - start, 2)} seconds ************************')

    return

if __name__ == '__main__':
    main()
