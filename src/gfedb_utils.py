import os
import sys
import logging
import re
import ast
import time
from Bio import AlignIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from seqann.models.annotation import Annotation
from Bio import SeqIO
from pyard import ARD
from csv import DictWriter
from pathlib import Path
from constants import *
from pympler import tracker, muppy, summary

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

tr = tracker.SummaryTracker()

expre_chars = ['N', 'Q', 'L', 'S']
isutr = lambda f: True if re.search("UTR", f) else False
to_second = lambda a: ":".join(a.split(":")[0:2]) + list(a)[-1] if list(a)[-1] in expre_chars and len(
    a.split(":")) > 2 else ":".join(a.split(":")[0:2])


def hla_alignments(dbversion):
    gen_aln = {l: {} for l in hla_loci}
    nuc_aln = {l: {} for l in hla_loci}
    prot_aln = {l: {} for l in hla_loci}

    #logging.info(f'HLA alignments:\n{hla_align}')

    for loc in hla_align:
        msf_gen = ''.join([data_dir, dbversion, "/", loc.split("-")[1], "_gen.msf"])
        msf_nuc = ''.join([data_dir, dbversion, "/", loc.split("-")[1], "_nuc.msf"])
        msf_prot = ''.join([data_dir, dbversion, "/", loc.split("-")[1], "_prot.msf"])

        logging.info(f'Loading {"/".join(msf_gen.split("/")[-3:])}')
        align_gen = AlignIO.read(open(msf_gen), "msf")
        gen_seq = {"HLA-" + a.name: str(a.seq) for a in align_gen}
        del align_gen
        logging.info(f'{str(len(gen_seq))} genomic alignments loaded')
        gen_aln.update({loc: gen_seq})

        logging.info(f'Loading {"/".join(msf_nuc.split("/")[-3:])}')
        align_nuc = AlignIO.read(open(msf_nuc), "msf")
        nuc_seq = {"HLA-" + a.name: str(a.seq) for a in align_nuc}
        del align_nuc
        logging.info(f'{str(len(nuc_seq))} nucelotide alignments loaded')
        nuc_aln.update({loc: nuc_seq})

        # https://github.com/ANHIG/IMGTHLA/issues/158
        # if str(dbversion) == ["3320", "3360"]:
        #    continue

        logging.info(f'Loading {"/".join(msf_prot.split("/")[-3:])}')
        align_prot = AlignIO.read(open(msf_prot), "msf")
        prot_seq = {"HLA-" + a.name: str(a.seq) for a in align_prot}
        del align_prot
        logging.info(f'{str(len(prot_seq))} protein alignments loaded')
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


# Streams dictionaries as rows to a file
def append_dict_as_row(file_path, dict_row):

    header = list(dict_row.keys())

    # Check if file exists
    csv_file = Path(file_path)
    if not csv_file.is_file():

        # Create the file and add the header
        with open(file_path, 'a+', newline='') as write_obj:
            dict_writer = DictWriter(write_obj, fieldnames=header)
            dict_writer.writeheader()

    # Do not add an else statement or the first line will be skipped
    with open(file_path, 'a+', newline='') as write_obj:
        dict_writer = DictWriter(write_obj, fieldnames=header)
        dict_writer.writerow(dict_row)

    return

# Outputs memory of objects during execution to sheck for memory leaks
def memory_profiler(mode='all'):

    # Print a summary of memory usage every n alleles
    all_objects = muppy.get_objects()
    sum2 = summary.summarize(all_objects)

    original_stdout = sys.stdout

    if mode == 'all' or mode == 'agg':
        with open("summary_agg.txt", "a+") as f:
            sys.stdout = f
            # tr.print_diff()
            summary.print_(sum2)
            sys.stdout = original_stdout;

    if mode == 'all' or mode == 'diff':
        with open("summary_diff.txt", "a+") as f:
            sys.stdout = f
            tr.print_diff()
            #summary.print_(sum2)
            sys.stdout = original_stdout;    

    return

# Build the datasets for the HLA graph
def build_hla_graph(**kwargs):

    dbversion, alignments, verbose, debug, gfe_maker, limit = \
        kwargs.get("dbversion"), \
        kwargs.get("alignments", False), \
        kwargs.get("verbose", False), \
        kwargs.get("debug", False), \
        kwargs.get("gfe_maker"), \
        kwargs.get("limit", None), \
    
    num_alleles = limit if limit else kwargs.get("num_alleles")

    def _stream_to_csv(a_gen, alignments, limit):

        i = 0
        total_time_elapsed = 0

        for idx, allele in enumerate(a_gen):

            start_time = time.time()

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

                    try:
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

                                gen_alignment = {
                                    "label": "GEN_ALIGN",
                                    "gfe_name": gfe,
                                    "hla_name": hla_name,
                                    "a_name": a_name, # hla_name.split("-")[1]
                                    "length": len(aligned_gen), # trim whitespace?
                                    "rank": "0", # TO DO: confirm how this value is derived
                                    "bp_sequence": aligned_gen,
                                    "aa_sequence": "",
                                    "imgt_release": imgt_release
                                }
                                
                                logging.info(f'Streaming genomic alignments to file...')
                                file_path = f'{data_dir}csv/all_alignments.{dbversion}.csv'
                                append_dict_as_row(file_path, gen_alignment)

                            if allele.description.split(",")[0] in nuc_aln[loc]:
                                aligned_nuc = nuc_aln[loc][allele.description.split(",")[
                                    0]]

                                nuc_alignment = {
                                    "label": "NUC_ALIGN",
                                    "gfe_name": gfe,
                                    "hla_name": hla_name,
                                    "a_name": a_name, # hla_name.split("-")[1]
                                    "length": len(aligned_nuc),
                                    "rank": "0", # TO DO: confirm how this value is derived
                                    "bp_sequence": aligned_nuc,
                                    "aa_sequence": "",
                                    "imgt_release": imgt_release
                                }

                                logging.info(f'Streaming nucleotide alignments to file...')
                                file_path = f'{data_dir}csv/all_alignments.{dbversion}.csv'
                                append_dict_as_row(file_path, nuc_alignment)

                            if allele.description.split(",")[0] in prot_aln[loc]:
                                aligned_prot = prot_aln[loc][allele.description.split(",")[
                                    0]]

                                prot_alignment = {
                                    "label": "PROT_ALIGN",
                                    "gfe_name": gfe,
                                    "hla_name": hla_name,
                                    "a_name": a_name, # hla_name.split("-")[1]
                                    "length": len(aligned_prot),
                                    "rank": "0", # TO DO: confirm how this value is derived
                                    "bp_sequence": "",
                                    "aa_sequence": aligned_prot,
                                    "imgt_release": imgt_release
                                }

                                logging.info(f'Streaming protein alignments to file...')
                                file_path = f'{data_dir}csv/all_alignments.{dbversion}.csv'
                                append_dict_as_row(file_path, prot_alignment)

                    except Exception as err:
                        logging.error(f'Failed to get data for allele ID {allele.id}')
                        logging.error(err)                        
                
                try:
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

                    logging.info(f'Streaming GFEs to file...')
                    file_name = ''.join([data_dir, "csv/gfe_sequence.{dbversion}.csv"])
                    append_dict_as_row(file_name, gfe_sequence)

                except Exception as err:
                    logging.error(f'Failed to write GFE data for allele ID {allele.id}')
                    logging.error(err)   

                try:
                    logging.info(f'Streaming groups to file...')
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

                        file_path = f'{data_dir}csv/all_groups.{dbversion}.csv'
                        append_dict_as_row(file_path, group_dict)

                    del groups

                except Exception as err:
                    logger.error(f'Failed to write groups for allele {allele.id}')
                    logger.error(err)

                try:
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

                    logging.info(f'Streaming CDS to file...')
                    file_path = f'{data_dir}csv/all_cds.{dbversion}.csv'
                    append_dict_as_row(file_path, cds)

                    del bp_seq
                    del aa_seq

                except Exception as err:
                    logger.error(f'Failed to write CDS data for allele {allele.id}')
                    logger.error(err)

                try:
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
                    
                    logging.info(f'Streaming features to file...')
                    for feature in features:
                        feature["gfe_name"] = gfe
                        feature["term"] = feature["term"].upper()
                        feature["allele_id"] = allele.id 
                        feature["hla_name"] = hla_name
                        feature["imgt_release"] = imgt_release

                        # Avoid null values in CSV for Neo4j import
                        feature["hash_code"] = "none" if not feature["hash_code"] else feature["hash_code"]

                        file_path = f'{data_dir}csv/all_features.{dbversion}.csv'
                        append_dict_as_row(file_path, feature)

                    del features

                except Exception as err:
                    logger.error(f'Failed to write features for allele {allele.id}')
                    logger.error(err)

            elapsed_time = time.time() - start_time
            alleles_remaining = num_alleles - (idx + 1)
            total_time_elapsed += elapsed_time
            avg_time_elapsed = total_time_elapsed / (idx + 1)
            #total_time_elapsed += ((alleles_remaining * elapsed_time) / 60)
            #avg_time_elapsed = total_time_elapsed / num_alleles
            #time_remaining = elapsed_time * alleles_remaining
            
            logging.info(f'Alleles processed: {idx + 1}')
            logging.info(f'Alleles remaining: {alleles_remaining}')
            logging.info(f'Elapsed time: {round(elapsed_time, 2)}')
            logging.info(f'Avg elapsed time: {round(avg_time_elapsed, 2)}')
            #logging.info(f'Estimated time remaining: {time.strftime("%H:%M:%S", time.gmtime(time_remaining))} minutes')
            
            # Break point for testing
            if limit and idx + 1 == limit:
                    break

            # Output memory profile to check for leaks
            if idx % 20 == 0:
                memory_profiler()

        return


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

    ### TO DO: move to build.sh
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

    logging.info("Streaming rows CSV files...")
    _stream_to_csv(
        a_gen=a_gen, 
        alignments=alignments, 
        limit=limit)

    return