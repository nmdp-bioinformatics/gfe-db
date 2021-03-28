#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TODO: **** ADD HAS_FEATURE
#       between SEQUENCE and features
import pandas as pd
from seqann.models.annotation import Annotation
from seqann.gfe import GFE
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pyard import ARD
from Bio import AlignIO
import logging
import argparse
import os
import sys
import urllib.request
import re
import ast
import time
#import pdb;
#from memory_profiler import profile
from pympler import tracker, muppy, summary
import gc

tr = tracker.SummaryTracker()

imgt_hla = 'https://www.ebi.ac.uk/ipd/imgt/hla/docs/release.html'
imgt_hla_media_url = 'https://media.githubusercontent.com/media/ANHIG/IMGTHLA/'
imgt_hla_raw_url = 'https://raw.githubusercontent.com/ANHIG/IMGTHLA/'

imgt_kir = 'https://www.ebi.ac.uk/ipd/kir/docs/version.html'
kir_url = 'ftp://ftp.ebi.ac.uk/pub/databases/ipd/kir/KIR.dat'

data_dir = os.path.dirname(__file__) + "/../data/"

logging.basicConfig(format='%(asctime)s - %(name)-25s - %(levelname)-5s - %(funcName)s:%(lineno)d: - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

expre_chars = ['N', 'Q', 'L', 'S']
isutr = lambda f: True if re.search("UTR", f) else False
to_second = lambda a: ":".join(a.split(":")[0:2]) + list(a)[-1] if list(a)[-1] in expre_chars and len(
    a.split(":")) > 2 else ":".join(a.split(":")[0:2])

# The alleles are removed when the allele_nodes.csv is built
skip_alleles = ["HLA-DRB5*01:11", "HLA-DRB5*01:12", "HLA-DRB5*01:13",
                "HLA-DRB5*02:03", "HLA-DRB5*02:04", "HLA-DRB5*02:05",
                "HLA-DRB5*01:01:02", "HLA-DRB5*01:03", "HLA-DRB5*01:05",
                "HLA-DRB5*01:06", "HLA-DRB5*01:07", "HLA-DRB5*01:09",
                "HLA-DRB5*01:10N", "HLA-C*05:208N", "HLA-C*05:206"]

hla_loci = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DQB1',
            'HLA-DPB1', 'HLA-DQA1', 'HLA-DPA1', 'HLA-DRB3',
            'HLA-DRB4', 'HLA-DRB5']

hla_align = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DQB1',
             'HLA-DPB1', 'HLA-DQA1', 'HLA-DPA1']

kir_loci = ["KIR3DS1", "KIR3DP1", "KIR3DL3", "KIR3DL2", "KIR3DL1",
            "KIR2DS5", "KIR2DS4", "KIR2DS3", "KIR2DS2", "KIR2DS1",
            "KIR2DP1", "KIR2DL5B", "KIR2DL5A", "KIR2DL4"]

kir_aligloci = ["KIR2DL4", "KIR2DP1", "KIR2DS1", "KIR2DS2", "KIR2DS3",
                "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3",
                "KIR3DP1"]

ard_groups = ['G', 'lg', 'lgx']

# # Not used
# def kir_alignments():
#     alignments = {l: {} for l in kir_loci}
#     for loc in kir_aligloci:
#         msf_file = data_dir + "/kir/" + loc + "_gen.msf"
#         align = AlignIO.read(open(msf_file), "msf")
#         names = {a.name: str(a.seq) for a in align}
#         alignments.update({loc: names})
#     return alignments


def hla_alignments(dbversion):
    gen_aln = {l: {} for l in hla_loci}
    nuc_aln = {l: {} for l in hla_loci}
    prot_aln = {l: {} for l in hla_loci}

    logging.info(f'HLA alignments:\n{hla_align}')

    for loc in hla_align:
        msf_gen = ''.join([data_dir, dbversion, "/", loc.split("-")[1], "_gen.msf"])
        msf_nuc = ''.join([data_dir, dbversion, "/", loc.split("-")[1], "_nuc.msf"])
        msf_prot = ''.join([data_dir, dbversion, "/", loc.split("-")[1], "_prot.msf"])

        logging.info(' '.join(["Loading", msf_gen]))
        align_gen = AlignIO.read(open(msf_gen), "msf")
        gen_seq = {"HLA-" + a.name: str(a.seq) for a in align_gen}
        # del align_gen.clear()
        logging.info(' '.join(["Loaded", str(len(gen_seq)), "genomic", loc, "alignments"]))
        gen_aln.update({loc: gen_seq})

        logging.info(' '.join(["Loading", msf_nuc]))
        align_nuc = AlignIO.read(open(msf_nuc), "msf")
        nuc_seq = {"HLA-" + a.name: str(a.seq) for a in align_nuc}
        # del align_nuc.clear()
        logging.info(' '.join(["Loaded", str(len(nuc_seq)), "nucleotide", loc, "alignments"]))
        nuc_aln.update({loc: nuc_seq})

        # https://github.com/ANHIG/IMGTHLA/issues/158
        # if str(dbversion) == ["3320", "3360"]:
        #    continue

        logging.info(' '.join(["Loading ", msf_prot]))
        align_prot = AlignIO.read(open(msf_prot), "msf")
        prot_seq = {"HLA-" + a.name: str(a.seq) for a in align_prot}
        # del align_prot.clear()
        logging.info(' '.join(["Loaded", str(len(prot_seq)), "protein", loc, "alignments"]))
        prot_aln.update({loc: prot_seq})

    return gen_aln, nuc_aln, prot_aln


def get_features(seqrecord):
    j = 3 if len(seqrecord.features) > 3 else len(seqrecord.features)
    fiveutr = [["five_prime_UTR", SeqRecord(seq=seqrecord.features[i].extract(seqrecord.seq), id="1")] for i in
               range(0, j) if seqrecord.features[i].type != "source"
               and seqrecord.features[i].type != "CDS" and isinstance(seqrecord.features[i], SeqFeature)
               and not seqrecord.features[i].qualifiers]
    feats = [[''.join([str(feat.type), "_", str(feat.qualifiers['number'][0])]), SeqRecord(seq=feat.extract(seqrecord.seq), id="1")]
             for feat in seqrecord.features if feat.type != "source"
             and feat.type != "CDS" and isinstance(feat, SeqFeature)
             and 'number' in feat.qualifiers]

    threeutr = []
    if len(seqrecord.features) > 1:
        threeutr = [["three_prime_UTR", SeqRecord(seq=seqrecord.features[i].extract(seqrecord.seq), id="1")] for i in
                    range(len(seqrecord.features) - 1, len(seqrecord.features)) if
                    seqrecord.features[i].type != "source"
                    and seqrecord.features[i].type != "CDS" and isinstance(seqrecord.features[i], SeqFeature)
                    and not seqrecord.features[i].qualifiers]

    feat_list = fiveutr + feats + threeutr
    annotation = {k[0]: k[1] for k in feat_list}

    return annotation

# Returns base pair and amino acid sequences from CDS data
def get_cds(allele):

    feat_types = [f.type for f in allele.features]
    bp_seq = None
    aa_seq = None
    
    if "CDS" in feat_types:
        cds_features = allele.features[feat_types.index("CDS")]
        if 'translation' in cds_features.qualifiers:

            if cds_features.location is None:
                logging.info(f"No CDS location for feature in allele: {allele.name}")
            else:
                bp_seq = str(cds_features.extract(allele.seq))
                aa_seq = cds_features.qualifiers['translation'][0]
                
    return bp_seq, aa_seq


# Build the datasets for the HLA graph
#@profile
def build_hla_graph(**kwargs):

    dbversion, alignments, verbose, debug, gfe_maker, limit, num_alleles = \
        kwargs.get("dbversion"), \
        kwargs.get("alignments", False), \
        kwargs.get("verbose", False), \
        kwargs.get("debug", False), \
        kwargs.get("gfe_maker"), \
        kwargs.get("limit", None), \
        kwargs.get("num_alleles")

    def _build_csv_files(a_gen, alignments, limit):

        i = 0

        ### Initialize lists for CSV output (input to LOAD CSV in Neo4j)
        # Lists contain unique dicts and are converted to dataframes, 
        # then output to CSV for Neo4j import
        gfe_sequences, gen_alignments, nuc_alignments, prot_alignments, all_features, all_groups, all_cds = \
            [], [], [], [], [], [], []


        #################
        logging.info(f'a_gen type: {type(a_gen)}')

        for idx, allele in enumerate(a_gen):

            t0 = time.time()
            # all_objects = muppy.get_objects()
            # sum1 = summary.summarize(all_objects)

            if hasattr(allele, 'seq'):
                hla_name = allele.description.split(",")[0]
                loc = allele.description.split(",")[0].split("*")[0]

                # if hla_name in skip_alleles:
                #     logging.info(
                #         "SKIPPING = " + allele.description.split(",")[0] + " " + dbversion)
                #     continue

                if debug and (loc != "HLA-A" and i > 20):
                    continue

                if (loc in hla_loci or loc == "DRB5") and (len(str(allele.seq)) > 5):
                    if debug:
                        logging.info(
                            "HLA = " + allele.description.split(",")[0] + " " + dbversion)

                    a_name = allele.description.split(",")[0].split("-")[1]
                    groups = [["HLA-" + ard.redux(a_name, grp), grp] if ard.redux(a_name, grp) != a_name else None for
                                grp in ard_groups]
                    seco = [[to_second(a_name), "2nd_FIELD"]]
                    groups = list(filter(None, groups)) + seco
                    complete_annotation = get_features(allele)
                    ann = Annotation(annotation=complete_annotation,
                                        method='match',
                                        complete_annotation=True)

                    # This process takes a long time
                    logging.info(f"Getting GFE data for allele {allele.id}...")
                    features, gfe = gfe_maker.get_gfe(ann, loc)

                    # gen_aln, nuc_aln, prot_aln
                    alignments_data = None
                    aligned_gen = ''
                    aligned_nuc = ''
                    aligned_prot = ''

                    if alignments:
                        if allele.description.split(",")[0] in gen_aln[loc]:
                            aligned_gen = gen_aln[loc][allele.description.split(",")[
                                0]]

                            # Separate CSV file, GFE foreign key: a_name
                            gen_alignment = {
                                "label": "GEN_ALIGN",
                                "gfe_name": gfe,
                                "hla_name": hla_name,
                                "a_name": a_name, # hla_name.split("-")[1]
                                "length": len(aligned_gen), # trim whitespace?
                                "rank": "0", # TO DO: confirm how this value is derived
                                "bp_sequence": aligned_gen,
                                "imgt_release": imgt_release
                            }                                    

                        if allele.description.split(",")[0] in nuc_aln[loc]:
                            aligned_nuc = nuc_aln[loc][allele.description.split(",")[
                                0]]

                            # Separate CSV file, GFE foreign key: a_name
                            nuc_alignment = {
                                "label": "NUC_ALIGN",
                                "gfe_name": gfe,
                                "hla_name": hla_name,
                                "a_name": a_name, # hla_name.split("-")[1]
                                "length": len(aligned_nuc),
                                "rank": "0", # TO DO: confirm how this value is derived
                                "bp_sequence": aligned_nuc,
                                "imgt_release": imgt_release
                            }

                        if allele.description.split(",")[0] in prot_aln[loc]:
                            aligned_prot = prot_aln[loc][allele.description.split(",")[
                                0]]

                            # Separate CSV file, GFE foreign key: a_name
                            prot_alignment = {
                                "label": "PROT_ALIGN",
                                "gfe_name": gfe,
                                "hla_name": hla_name,
                                "a_name": a_name, # hla_name.split("-")[1]
                                "length": len(aligned_prot),
                                "rank": "0", # TO DO: confirm how this value is derived
                                "aa_sequence": aligned_prot,
                                "imgt_release": imgt_release
                            }

                        # Append data to respective list
                        alignments_data = zip(
                            [gen_alignments, nuc_alignments, prot_alignments],
                            [gen_alignment, nuc_alignment, prot_alignment]
                        )

                        for _list, _dict in alignments_data:
                            _list.append(_dict)
                        
               

                ### Build dicts describing nodes and edges for each allele
                # Separate CSV file
                gfe_sequence = {
                    "gfe_name": gfe,
                    "allele_id": allele.id,
                    "locus": loc,
                    "hla_name": hla_name,
                    "a_name": a_name, # hla_name.split("-")[1]
                    "sequence": str(allele.seq),
                    "length": len(str(allele.seq)),
                    "imgt_release": imgt_release
                }

                # Separate CSV file, GFE foreign key: allele_id
                allele_groups = []

                for group in groups:
                    group_dict = {
                        "gfe_name": gfe,
                        "allele_id": allele.id,
                        "hla_name": hla_name,
                        "a_name": a_name,
                        "ard_id": group[0],
                        "ard_name": group[1],
                        "locus": loc,
                        "imgt_release": imgt_release
                    }

                    allele_groups.append(group_dict)

                # Build CDS dict for CSV export, foreign key: allele_id, hla_name
                bp_seq, aa_seq = get_cds(allele)

                cds = {
                    "gfe_name": gfe,
                    "allele_id": allele.id,
                    "hla_name": hla_name,
                    "bp_sequence": bp_seq,
                    "aa_sequence": aa_seq,
                    "imgt_release": imgt_release
                }

                # features preprocessing steps
                # 1) Convert seqann type to python dict using literal_eval
                # 2) add GFE foreign keys: allele_id, hla_name
                # 3) calculate columns: length

                # features contains list of seqann objects, converts to dict, destructive step
                features = \
                    [ast.literal_eval(str(feature) \
                        .replace('\'', '"') \
                        .replace('\n', '')) \
                        for feature in features]               

                # Append allele id's
                # Note: Some alleles may have the same feature, but it may not be the same rank, 
                # so a feature should be identified with its allele by allele_id or HLA name
                for feature in features:
                    feature["gfe_name"] = gfe
                    feature["term"] = feature["term"].upper()
                    feature["allele_id"] = allele.id 
                    feature["hla_name"] = hla_name
                    feature["imgt_release"] = imgt_release

                    # Avoid null values in CSV for Neo4j import
                    feature["hash_code"] = "none" if not feature["hash_code"] else feature["hash_code"]

                # Append data to respective list
                data = zip(
                    [gfe_sequences, all_cds],
                    [gfe_sequence, cds]
                )

                for _list, _dict in data:
                    _list.append(_dict)

                # Alignments, features, and ARD groups can all be concatenated since the keys are the same
                if alignments_data:
                    all_alignments = gen_alignments + nuc_alignments + prot_alignments
                    # del alignments_data

                    # Drop duplicate rows
                    all_alignments = pd.DataFrame(all_alignments) \
                        .drop_duplicates() \
                        .to_dict('records')

                all_features = all_features + features        
                all_groups = all_groups + allele_groups

            elapsed_time = round(time.time()-t0, 2)
            time_remaining = round(((num_alleles - idx) * elapsed_time) / 60, 2)

            alleles_remaining = num_alleles - idx

            logging.info(f'Alleles processed: {idx}')
            logging.info(f'Alleles remaining: {alleles_remaining}')
            logging.info(f'Elapsed time: {elapsed_time}')
            logging.info(f'Estimated time remaining: {time_remaining} minutes')
            
            # Break point for testing
            if limit and idx == limit:
                    break

            # Print a summary of memory usage every n alleles
            if idx % 20 == 0:
                all_objects = muppy.get_objects()
                sum2 = summary.summarize(all_objects)
                # #diff = summary.get_diff(sum1, sum2)
                # #summary.print_(diff)
                # #logging.info(f'Memory usage:\n{}\n')
                # summary.print_(sum2)
                # with open("summary.json", "a+") as f:
                #     f.write(json.dumps(sum1))

                original_stdout = sys.stdout
                with open("summary_agg.txt", "a+") as f:
                    sys.stdout = f
                    # tr.print_diff()
                    summary.print_(sum2)
                    sys.stdout = original_stdout;

                with open("summary_diff.txt", "a+") as f:
                    sys.stdout = f
                    tr.print_diff()
                    #summary.print_(sum2)
                    sys.stdout = original_stdout;                
                
        
        #################

        
        csv_output = {
            "gfe_sequences": gfe_sequences,
            "all_features": all_features,
            "all_groups": all_groups,
            "all_cds": all_cds
        }

        # Add alignments data if there is
        if alignments:
            csv_output["all_alignments"] = all_alignments
            # del all_alignments

        # # To Do : use .clear() method for dictionaries
        # del allele_groups
        # del features
        # del gfe_sequence
        # del cds
        # del gfe_sequences
        # del all_features
        # del all_groups
        # del all_cds

        gc.collect()

        return csv_output


    # Loop through DB versions and build CSVs
    #for dbversion in dbversions:

    imgt_release = f'{dbversion[0]}.{dbversion[1:3]}.{dbversion[3]}'
    db_striped = ''.join(dbversion.split("."))
    
    logging.debug(f'dbversion: {dbversion}')
    logging.debug(f'imgt_release: {imgt_release}')
    logging.debug(f'db_striped: {db_striped}')

    if alignments:
        gen_aln, nuc_aln, prot_aln = hla_alignments(db_striped)

    logging.info("Loading ARD...")
    ard = ARD(db_striped)

    # The github URL changed from 3350 to media
    if int(db_striped) < 3350:
        dat_url = ''.join([imgt_hla_raw_url, db_striped, '/hla.dat'])
    else:
        dat_url = ''.join([imgt_hla_media_url, db_striped, '/hla.dat'])

    dat_file = ''.join([data_dir, 'hla.', db_striped, ".dat"])

    # Downloading DAT file
    logging.info("Downloading DAT file...")
    if not os.path.isfile(dat_file):
        if verbose:
            logging.info("Downloading dat file from " + dat_url)
        urllib.request.urlretrieve(dat_url, dat_file)

    # Parse DAT file
    logging.info("Parsing DAT file...")
    a_gen = SeqIO.parse(dat_file, "imgt")

    if verbose:
        logging.info("Finished parsing dat file")

    
    logging.info("Building CSV files...")
    csv_output = \
        _build_csv_files(
            a_gen=a_gen, 
            alignments=alignments, 
            limit=limit)

    # del a_gen
    # del alignments

    return csv_output


# Write CSV files to local directory
def write_csv_file(csv_output, dbversion, path):
    """Takes a dict of form "csv_name": data where csv_name is the CSV file to export
    and data is a list of dictionaries.
    """

    try:
        # Output to CSV, include dbversion in name
        for csv_name, data in csv_output.items():
            dataframe = pd.DataFrame(data)
            file_name = path + f"{csv_name}.{dbversion}.csv"
            dataframe.to_csv(file_name, index=False)
            logging.info(f'Saved {file_name.split("/")[-1]} to {path}')

        return 

    except Exception as err:
        logging.error(f'Failed to save CSV files to "{path}"')
        raise err


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

    parser.add_argument("-r", "--releases",
                        required=False,
                        help="IMGT/DB releases",
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

    logging.info(f'args:\n{vars(args)}')

    # Not used
    # out_dir = args.out_dir if args.out_dir else ""

    release_n = args.number
    releases = args.releases if args.releases else None
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

    # Get last five IMGT/HLA releases
    if releases:
        dbversions = [db for db in releases.split(",")]
    else:
        dbversions = pd.read_html(imgt_hla)[0]['Release'][0:release_n].tolist()

    # # Get latest IMGT/KIR release
    # kir_release = pd.read_html(imgt_kir)[0]['Release'][0]

    gfe_maker = GFE(verbose=verbose, 
        verbosity=verbosity,
        load_features=False, 
        store_features=True,
        loci=load_loci)

    if kir:

        # # Replace this block with function:
        # build_kir_graph(**kwargs)

        logging.info(f'KIR: {kir}')

        # if verbose:
        #     logging.info("Adding KIR to GFE DB")

        # kir_file = data_dir + 'KIR.dat'

        # if alignments:
        #     aligned = kir_alignments()

        # # Downloading KIR
        # if not os.path.isfile(kir_file):
        #     if verbose:
        #         logging.info("Downloading KIR dat file from " + kir_url)
        #     urllib.request.urlretrieve(kir_url, kir_file)

        # kir_gen = SeqIO.parse(kir_file, "imgt")
        # if verbose:
        #     logging.info("Finished parsing KIR dat file")

        # i = 0
        # for allele in kir_gen:
        #     if hasattr(allele, 'seq'):
        #         loc = allele.description.split(",")[0].split("*")[0]
        #         if loc in kir_loci and len(str(allele.seq)) > 5:
        #             if debug:
        #                 logging.info("KIR = " + allele.description.split(",")[0] + " " + kir_release)

        #             groups = []
        #             complete_annotation = get_features(allele)
        #             ambigs = [a for a in complete_annotation if re.search("/", a)]

        #             aligned_seq = ''
        #             if alignments:
        #                 if allele.description.split(",")[0] in aligned[loc]:
        #                     aligned_seq = aligned[loc][allele.description.split(",")[0]]

        #             if ambigs:
        #                 logging.info("AMBIGS " + allele.description.split(",")[0] + " " + kir_release)
        #                 annotations = []
        #                 for ambig in ambigs:
        #                     if debug:
        #                         logging.info("AMBIG = " + ambig)
        #                     aterm = ambig.split("/")[0].split("_")[0]
        #                     anno = {a: complete_annotation[a] for a in complete_annotation if a not in ambigs}
        #                     anno.update({ambig.split("/")[0]: complete_annotation[ambig]})
        #                     annotations.append(anno)

        #                     anno2 = {a: complete_annotation[a] for a in complete_annotation if a not in ambigs}
        #                     anno2.update({aterm + "_" + ambig.split("/")[1]: complete_annotation[ambig]})
        #                     annotations.append(anno2)

        #                 for annotation in annotations:
        #                     ann = Annotation(annotation=annotation,
        #                                      method='match',
        #                                      complete_annotation=True)

        #                     features, gfe = gfe_maker.get_gfe(ann, loc)
        #                     
        #                     # Build the KIR graph
        #                     build_kir_graph(...)

        #             else:
        #                 ann = Annotation(annotation=complete_annotation,
        #                                  method='match',
        #                                  complete_annotation=True)
        #                 features, gfe = gfe_maker.get_gfe(ann, loc)

        #                 # Build the KIR graph
        #                 build_kir_graph(...)

    # TO DO: for this for loop into the build.sh script
    for dbversion in dbversions:
        
        logging.info(f'\n\nBuilding graph for IMGTHLA version {dbversion[0]}.{dbversion[1:3]}.{dbversion[3]}...')
        logging.info(f'Limit: {args.limit}')
        csv_output = build_hla_graph(
            dbversion=dbversion, 
            alignments=align, 
            verbose=verbose,
            limit=args.limit,
            num_alleles=num_alleles,
            gfe_maker=gfe_maker)

        write_csv_file(csv_output, dbversion, path=''.join([data_dir, "csv/"]))

        print("csv_output: ", type(csv_output))
        print("csv_output keys: ", csv_output.keys())
        logging.info(f'Finished build for version {dbversion[0]}.{dbversion[1:3]}.{dbversion[3]}')

    logging.info(f'****** Build complete ******')

    return

if __name__ == '__main__':
    """The following will be run if file is executed directly,
    but not if imported as a module"""
    main()
