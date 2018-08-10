#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TODO: **** ADD HAS_FEATURE
#       between SEQUENCE and features
import pandas as pd
from seqann.models.annotation import Annotation
from Bio.SeqFeature import SeqFeature
from pygfe.pygfe import pyGFE
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pyard import ARD
from Bio import AlignIO
import logging
import argparse
import sys
import os
import urllib.request
import re


imgt_hla = "https://www.ebi.ac.uk/ipd/imgt/hla/docs/release.html"
imgt_kir = "https://www.ebi.ac.uk/ipd/kir/docs/version.html"
kir_url = "ftp://ftp.ebi.ac.uk/pub/databases/ipd/kir/KIR.dat"

logging.basicConfig(format='%(asctime)s - %(name)-35s - %(levelname)-5s - %(funcName)s %(lineno)d: - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    level=logging.INFO)

expre_chars = ['N', 'Q', 'L', 'S']
isutr = lambda f: True if re.search("UTR", f) else False
to_second = lambda a: ":".join(a.split(":")[0:2]) + list(a)[-1] if list(a)[-1] in expre_chars and len(a.split(":")) > 2 else ":".join(a.split(":")[0:2])

lastseqid = 1
lastid = 1
lastcdsid = 1

seqids = {}
cdsids = {}
alleleids = {}
group_edges = {}
trans_edges = {}

skip_alleles = ["HLA-DRB5*01:11", "HLA-DRB5*01:12", "HLA-DRB5*01:13",
                "HLA-DRB5*02:03", "HLA-DRB5*02:04", "HLA-DRB5*02:05",
                "HLA-DRB5*01:01:02", "HLA-DRB5*01:03", "HLA-DRB5*01:05",
                "HLA-DRB5*01:06", "HLA-DRB5*01:07", "HLA-DRB5*01:09",
                "HLA-DRB5*01:10N"]

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


def kir_alignments():

    alignments = {l: {} for l in kir_loci}
    data_dir = os.path.dirname(__file__)
    for loc in kir_aligloci:
        sth_file = data_dir + "/../data/" + loc + "_gen.sth"
        align = AlignIO.read(open(sth_file), "stockholm")
        names = {a.name: str(a.seq) for a in align}
        alignments.update({loc: names})
    return alignments


def hla_alignments(dbversion):

    gen_aln = {l: {} for l in hla_loci}
    nuc_aln = {l: {} for l in hla_loci}
    prot_aln = {l: {} for l in hla_loci}

    data_dir = os.path.dirname(__file__)
    for loc in hla_align:
        sth_gen = data_dir + "/../data/" + dbversion + "/" + loc.split("-")[1] + "_gen.sth"
        sth_nuc = data_dir + "/../data/" + dbversion + "/" + loc.split("-")[1] + "_nuc.sth"
        sth_prot = data_dir + "/../data/" + dbversion + "/" + loc.split("-")[1] + "_prot.sth"

        logging.info("Loading " + sth_gen)
        align_gen = AlignIO.read(open(sth_gen), "stockholm")
        gen_seq = {"HLA-" + a.name: str(a.seq) for a in align_gen}
        logging.info("Loaded " + str(len(gen_seq)) + " genomic " + loc + " sequences")
        gen_aln.update({loc: gen_seq})

        logging.info("Loading " + sth_nuc)
        align_nuc = AlignIO.read(open(sth_nuc), "stockholm")
        nuc_seq = {"HLA-" + a.name: str(a.seq) for a in align_nuc}
        logging.info("Loaded " + str(len(nuc_seq)) + " nuc " + loc + " sequences")
        nuc_aln.update({loc: nuc_seq})

        # https://github.com/ANHIG/IMGTHLA/issues/158 
        if str(dbversion) == "3320":
            continue

        logging.info("Loading " + sth_prot)
        align_prot = AlignIO.read(open(sth_prot), "stockholm")
        prot_seq = {"HLA-" + a.name: str(a.seq) for a in align_prot}
        logging.info("Loaded " + str(len(prot_seq)) + " prot " + loc + " sequences")
        prot_aln.update({loc: prot_seq})

    return gen_aln, nuc_aln, prot_aln


def get_features(seqrecord):
    j = 3 if len(seqrecord.features) > 3 else len(seqrecord.features)
    fiveutr = [["five_prime_UTR", SeqRecord(seq=seqrecord.features[i].extract(seqrecord.seq), id="1")] for i in range(0, j) if seqrecord.features[i].type != "source"
               and seqrecord.features[i].type != "CDS" and isinstance(seqrecord.features[i], SeqFeature)
               and not seqrecord.features[i].qualifiers]
    feats = [[str(feat.type + "_" + feat.qualifiers['number'][0]), SeqRecord(seq=feat.extract(seqrecord.seq), id="1")]
             for feat in seqrecord.features if feat.type != "source"
             and feat.type != "CDS" and isinstance(feat, SeqFeature)
             and 'number' in feat.qualifiers]

    threeutr = []
    if len(seqrecord.features) > 1:
        threeutr = [["three_prime_UTR", SeqRecord(seq=seqrecord.features[i].extract(seqrecord.seq), id="1")] for i in range(len(seqrecord.features)-1, len(seqrecord.features)) if seqrecord.features[i].type != "source"
                    and seqrecord.features[i].type != "CDS" and isinstance(seqrecord.features[i], SeqFeature)
                    and not seqrecord.features[i].qualifiers]

    feat_list = fiveutr + feats + threeutr
    annotation = {k[0]: k[1] for k in feat_list}

    return(annotation)


def build_graph(groups, gfe, allele, features, dbversion,
                align_gen, align_nuc, align_prot, allele_type, load_alignments):

    global lastseqid
    global lastid
    global lastcdsid

    loc = gfe.split("w")[0]
    full_seq = str(allele.seq)

    uniqseq = True
    uniqalignp = True
    uniqaligng = True
    uniqalignn = True
    uniqallele = True

    seq_edges = []
    seq_nodes = []
    grp_edges = []
    cds_nodes = []
    group_nodes = []

    if full_seq in seqids:
        uniqseq = False
        fullseqid = seqids[full_seq]
    else:
        fullseqid = lastseqid
        seqids.update({full_seq: lastseqid})
        lastseqid += 1

    if load_alignments:
        # Getting alignment IDs
        if not align_gen:
            uniqaligng = False
        else:
            if align_gen in seqids:
                uniqaligng = False
                genid = seqids[align_gen]
            else:
                genid = lastseqid
                seqids.update({align_gen: lastseqid})
                lastseqid += 1

        if not align_nuc:
            uniqalignn = False
        else:
            if align_nuc in seqids:
                uniqalignn = False
                nucid = seqids[align_nuc]
            else:
                nucid = lastseqid
                seqids.update({align_nuc: lastseqid})
                lastseqid += 1

        if not align_prot:
            uniqalignp = False
        else:
            if align_prot in seqids:
                uniqalignp = False
                protid = seqids[align_prot]
            else:
                protid = lastseqid
                seqids.update({align_prot: lastseqid})
                lastseqid += 1

    # GFE ID
    if gfe in alleleids:
        gfeid = alleleids[gfe]
        gfenode = []
    else:
        gfeid = lastid
        alleleids.update({gfe: gfeid})
        lastid += 1
        gfenode = [[gfeid, gfe, "GFE", loc]]

    # Allele ID
    allele_name = allele.description.split(",")[0]
    if allele_name in alleleids:
        uniqallele = False
        alleleid = alleleids[allele_name]
    else:
        alleleid = lastid
        alleleids.update({allele_name: alleleid})
        lastid += 1

    # GFE
    seq_edges.append([gfeid, fullseqid, dbversion, 0, "HAS_SEQUENCE"])

    # Alleles
    seq_edges.append([alleleid, fullseqid, dbversion, 0, "HAS_SEQUENCE"])

    if load_alignments:
        if align_gen:
            seq_edges.append([gfeid, genid, dbversion, 0, "HAS_ALIGNMENT"])
            seq_edges.append([alleleid, genid, dbversion, 0, "HAS_ALIGNMENT"])

        if align_nuc:
            seq_edges.append([gfeid, nucid, dbversion, 0, "HAS_ALIGNMENT"])
            seq_edges.append([alleleid, nucid, dbversion, 0, "HAS_ALIGNMENT"])

        if align_prot:
            seq_edges.append([gfeid, protid, dbversion, 0, "HAS_ALIGNMENT"])
            seq_edges.append([alleleid, protid, dbversion, 0, "HAS_ALIGNMENT"])

    if groups:
        for g in groups:
            grp = g[0]
            g_type = g[1]
            grp_n = "-".join([grp, g_type])
            if grp_n in alleleids:
                gid = alleleids[grp_n]
            else:
                gid = lastid
                alleleids.update({grp_n: gid})
                group_nodes.append([gid, grp, g_type, loc])
                lastid += 1

            allele_grp = "-".join([str(alleleid), str(gid)])
            if allele_grp not in group_edges:
                g_edge = [gid, alleleid, dbversion, g_type]
                grp_edges.append(g_edge)
                group_edges.update({allele_grp: g_edge})

    ## Only if GFE DOESN"T EXIST
    gfeedge = [[alleleid, gfeid, dbversion, "HAS_GFE"]]

    allelenode = gfenode
    if uniqallele:
        allelenode = allelenode + [[alleleid, allele.description.split(",")[0],
                                   allele_type, loc]]
        if group_nodes:
            allelenode = allelenode + group_nodes

    #alleleId:ID(ALLELE),name,alleletype:LABEL,locus,imgtdb
    #sequenceId:ID(SEQUENCE),sequence,name,feature:LABEL,rank,length,nuc:string[]

    feat_types = [f.type for f in allele.features]

    # *** NO CDS or translation
    # Have them be separate nodes
    tr_edges = []
    if "CDS" in feat_types:
        if 'translation' in allele.features[feat_types.index("CDS")].qualifiers:
            if uniqseq:
                seq_nodes.append([fullseqid, full_seq, "SEQUENCE", "SEQUENCE", 0,
                                  len(full_seq)])

            if load_alignments:
                if uniqaligng:
                    seq_nodes.append([genid, '', "GEN_ALIGN", "GEN_ALIGN", 0,
                                      len(align_gen), ";".join(list(align_gen))])

                if uniqalignn:
                    seq_nodes.append([nucid, '', "NUC_ALIGN", "NUC_ALIGN", 0,
                                      len(align_nuc), ";".join(list(align_nuc))])

                if uniqalignp:
                    seq_nodes.append([protid, '', "PROT_ALIGN", "PROT_ALIGN", 0,
                                      len(align_prot), ";".join(list(align_prot))])

            #sequenceId:ID(SEQUENCE),sequence,name,feature:LABEL,rank,length,nuc:string[]
            cds_seq = allele.features[feat_types.index("CDS")].extract(allele.seq)
            tran_seq = allele.features[feat_types.index("CDS")].qualifiers['translation'][0]

            if cds_seq in cdsids:
                cdsseqid = cdsids[cds_seq]
            else:
                cdsseqid = lastcdsid
                cdsids.update({cds_seq: lastcdsid})
                lastcdsid += 1
                cds_nodes.append([cdsseqid, "CDS", "CDS", cds_seq,  tran_seq])

            tr_edges = []
            seq_cds = "-".join([str(fullseqid), str(cdsseqid)])
            if seq_cds not in trans_edges:
                tr_edges = [[fullseqid, cdsseqid, "HAS_CDS"]]
                trans_edges.update({seq_cds: tr_edges})

        else:
            if uniqseq:
                seq_nodes.append([fullseqid, full_seq, "SEQUENCE", "SEQUENCE", 0,
                                  len(full_seq), ''])
            if load_alignments:
                if uniqaligng:
                    seq_nodes.append([genid, '', "GEN_ALIGN", "GEN_ALIGN", 0,
                                      len(align_gen), ";".join(list(align_gen))])

                if uniqalignn:
                    seq_nodes.append([nucid, '', "NUC_ALIGN", "NUC_ALIGN", 0,
                                      len(align_nuc), ";".join(list(align_nuc))])

                if uniqalignp:
                    seq_nodes.append([protid, '', "PROT_ALIGN", "PROT_ALIGN", 0,
                                      len(align_prot), ";".join(list(align_prot))])

            cds_seq = allele.features[feat_types.index("CDS")].extract(allele.seq)

            if cds_seq in cdsids:
                cdsseqid = cdsids[cds_seq]
            else:
                cdsseqid = lastcdsid
                cdsids.update({cds_seq: lastcdsid})
                lastcdsid += 1

                cds_nodes.append([cdsseqid, "CDS", "CDS", cds_seq, ''])

    else:
        if uniqseq:
            seq_nodes.append([fullseqid, full_seq, "SEQUENCE", "SEQUENCE", 0,
                              len(full_seq), ''])
        
        if load_alignments:
            if uniqaligng:
                seq_nodes.append([genid, '', "GEN_ALIGN", "GEN_ALIGN", 0,
                                  len(align_gen), ";".join(list(align_gen))])

            if uniqalignn:
                seq_nodes.append([nucid, '', "NUC_ALIGN", "NUC_ALIGN", 0,
                                  len(align_nuc), ";".join(list(align_nuc))])

            if uniqalignp:
                seq_nodes.append([protid, '', "PROT_ALIGN", "PROT_ALIGN", 0,
                                  len(align_prot), ";".join(list(align_prot))])

    for feat in features:
        feat_type = feat.term.upper()
        feat_seq = feat.sequence
        seq_type = "-".join([feat_seq, feat_type, str(feat.rank)])
        if seq_type in seqids:
            seqid = seqids[seq_type]
            seq_edges.append([alleleid, seqid, dbversion,
                              feat.accession, "HAS_FEATURE"])
            seq_edges.append([gfeid, seqid, dbversion,
                              feat.accession, "HAS_FEATURE"])
        else:
            seqid = lastseqid
            seqids.update({seq_type: lastseqid})
            lastseqid += 1

            seq_edges.append([alleleid, seqid, dbversion,
                              feat.accession, "HAS_FEATURE"])
            seq_edges.append([gfeid, seqid, dbversion,
                              feat.accession, "HAS_FEATURE"])
            if uniqseq:
                seq_nodes.append([seqid, feat_seq, feat_type, "FEATURE",
                                  feat.rank, len(feat_seq), ''])

    return allelenode, gfeedge, seq_nodes, cds_nodes, seq_edges, tr_edges, grp_edges


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
                        action='store_true')

    parser.add_argument("-d", "--debug",
                        required=False,
                        help="Bool for debugging",
                        action='store_true')

    parser.add_argument("-o", "--outdir",
                        required=True,
                        help="Output directory",
                        type=str)

    parser.add_argument("-n", "--number",
                        required=False,
                        help="Number of IMGT/DB releases",
                        default=1,
                        type=int)

    parser.add_argument("-r", "--releases",
                        required=False,
                        help="IMGT/DB releases",
                        type=str)

    parser.add_argument("-v", "--verbose",
                        help="Option for running in verbose",
                        action='store_true')

    data_dir = os.path.dirname(__file__)
    args = parser.parse_args()

    outdir = args.outdir

    load_loci = hla_loci + kir_loci
    release_n = args.number
    releases = args.releases
    verbosity = 1

    align = False
    kir = False
    debug = False
    verbose = False

    if args.kir:
        kir = True

    if args.align:
        align = True

    if args.verbose:
        verbose = True

    if kir:
        load_loci = hla_loci + kir_loci
    else:
        load_loci = hla_loci

    if args.debug:
        logging.info("Running in debug mode")
        load_loci = ["HLA-A"]
        kir = False
        debug = True
        verbose = True
        verbosity = 2
        release_n = 1

    gfe_e = []
    seq_e = []
    seq_n = []
    cds_n = []
    grp_e = []
    trs_e = []
    allele_n = []

    # Get last five IMGT/HLA releases
    if releases:
        dbversions = [db for db in releases.split(",")]
    else:
        dbversions = pd.read_html(imgt_hla)[0]['Release'][0:release_n].tolist()

    # Get lastest IMGT/KIR release
    kir_release = pd.read_html(imgt_kir)[0][0][1]

    gfe_maker = pyGFE(verbose=verbose, verbosity=verbosity,
                      load_features=True, store_features=True,
                      loci=load_loci)

    if kir:
        if verbose:
            logging.info("Adding KIR to GFE DB")

        kir_file = data_dir + '/../data/KIR.dat'

        if align:
            aligned = kir_alignments()

        # Downloading KIR
        if not os.path.isfile(kir_file):
            if verbose:
                logging.info("Downloading KIR dat file from " + kir_url)
            urllib.request.urlretrieve(kir_url, kir_file)

        kir_gen = SeqIO.parse(kir_file, "imgt")
        if verbose:
            logging.info("Finished parsing KIR dat file")

        i = 0
        for allele in kir_gen:
            if hasattr(allele, 'seq'):
                loc = allele.description.split(",")[0].split("*")[0]
                if loc in kir_loci and len(str(allele.seq)) > 5:
                    if verbose:
                        logging.info("KIR = " + allele.description.split(",")[0] + " " + kir_release)

                    groups = []
                    complete_annotation = get_features(allele)
                    ambigs = [a for a in complete_annotation if re.search("/", a)]

                    aligned_seq = ''
                    if align:
                        if allele.description.split(",")[0] in aligned[loc]:
                            aligned_seq = aligned[loc][allele.description.split(",")[0]]

                    if ambigs:
                        logging.info("AMBIGS " + allele.description.split(",")[0] + " " + kir_release)
                        annotations = []
                        for ambig in ambigs:
                            logging.info("AMBIG = " + ambig)
                            aterm = ambig.split("/")[0].split("_")[0]
                            anno = {a: complete_annotation[a] for a in complete_annotation if a not in ambigs}
                            anno.update({ambig.split("/")[0]: complete_annotation[ambig]})
                            annotations.append(anno)

                            anno2 = {a: complete_annotation[a] for a in complete_annotation if a not in ambigs}
                            anno2.update({aterm + "_" + ambig.split("/")[1]: complete_annotation[ambig]})
                            annotations.append(anno2)

                        for annotation in annotations:
                            ann = Annotation(annotation=annotation,
                                             method='match',
                                             complete_annotation=True)

                            features, gfe = gfe_maker.get_gfe(ann, loc)
                            (allelenode, gfeedge,
                             seq_nodes, cds_nodes, seq_edges,
                             trans_edge, grp_edges) = build_graph(groups,
                                                                  gfe,
                                                                  allele,
                                                                  features,
                                                                  kir_release,
                                                                  aligned_seq,
                                                                  '',
                                                                  '',
                                                                  "IMGT_KIR",
                                                                  align)

                            gfe_e += gfeedge
                            seq_e += seq_edges
                            seq_n += seq_nodes
                            allele_n += allelenode
                            grp_e += grp_edges
                            trs_e += trans_edge
                            cds_n += cds_nodes
                        i += 1

                    else:
                        ann = Annotation(annotation=complete_annotation,
                                         method='match',
                                         complete_annotation=True)
                        features, gfe = gfe_maker.get_gfe(ann, loc)

                        (allelenode, gfeedge,
                         seq_nodes, cds_nodes, seq_edges,
                         trans_edge, grp_edges) = build_graph(groups,
                                                              gfe,
                                                              allele,
                                                              features,
                                                              kir_release,
                                                              aligned_seq,
                                                              '',
                                                              '',
                                                              "IMGT_KIR",
                                                              align)

                        gfe_e += gfeedge
                        seq_e += seq_edges
                        seq_n += seq_nodes
                        allele_n += allelenode
                        grp_e += grp_edges
                        trs_e += trans_edge
                        cds_n += cds_nodes
                        i += 1

    # Loop through DB versions
    for dbversion in dbversions:

        db_striped = ''.join(dbversion.split("."))

        if align:
            gen_aln, nuc_aln, prot_aln = hla_alignments(db_striped)

        ard = ARD(db_striped)

        dat_url = 'https://raw.githubusercontent.com/ANHIG/IMGTHLA/' \
                  + db_striped + '/hla.dat'
        dat_file = data_dir + '/hla.' + str(db_striped) + ".dat"

        # Downloading DAT file
        if not os.path.isfile(dat_file):
            if verbose:
                logging.info("Downloading dat file from " + dat_url)
            urllib.request.urlretrieve(dat_url, dat_file)

        cmd = "perl -p -i -e 's/[^\\x00-\\x7F]//g' " + dat_file
        os.system(cmd)

        a_gen = SeqIO.parse(dat_file, "imgt")
        if verbose:
            logging.info("Finished parsing dat file")

        i = 0
        for allele in a_gen:
            if hasattr(allele, 'seq'):
                hla_name = allele.description.split(",")[0]
                loc = allele.description.split(",")[0].split("*")[0]
                if hla_name in skip_alleles:
                    logging.info("SKIPPING = " + allele.description.split(",")[0] + " " + dbversion)
                    continue

                if (debug and (loc != "HLA-A" and i > 20)):
                    continue

                if (loc in hla_loci or loc == "DRB5") and (len(str(allele.seq)) > 5):
                    if verbose:
                        logging.info("HLA = " + allele.description.split(",")[0] + " " + dbversion)

                    a_name = allele.description.split(",")[0].split("-")[1]
                    groups = [["HLA-" + ard.redux(a_name, grp), grp] if ard.redux(a_name, grp) != a_name else None for grp in ard_groups]
                    seco = [[to_second(a_name), "2nd_FIELD"]]
                    groups = list(filter(None, groups)) + seco
                    complete_annotation = get_features(allele)
                    ann = Annotation(annotation=complete_annotation,
                                     method='match',
                                     complete_annotation=True)
                    features, gfe = gfe_maker.get_gfe(ann, loc)

                    #gen_aln, nuc_aln, prot_aln
                    aligned_gen = ''
                    aligned_nuc = ''
                    aligned_prot = ''

                    if align:
                        if allele.description.split(",")[0] in gen_aln[loc]:
                            aligned_gen = gen_aln[loc][allele.description.split(",")[0]]

                        if allele.description.split(",")[0] in nuc_aln[loc]:
                            aligned_nuc = nuc_aln[loc][allele.description.split(",")[0]]

                        if allele.description.split(",")[0] in prot_aln[loc]:
                            aligned_prot = prot_aln[loc][allele.description.split(",")[0]]

                    (allelenode, gfeedge,
                     seq_nodes, cds_nodes, seq_edges,
                     trans_edge, grp_edges) = build_graph(groups,
                                                          gfe,
                                                          allele,
                                                          features,
                                                          dbversion,
                                                          aligned_gen,
                                                          aligned_nuc,
                                                          aligned_prot,
                                                          "IMGT_HLA",
                                                          align)

                    gfe_e += gfeedge
                    seq_e += seq_edges
                    seq_n += seq_nodes
                    allele_n += allelenode
                    grp_e += grp_edges
                    trs_e += trans_edge
                    cds_n += cds_nodes
                    i += 1
        if verbose:
            logging.info("Finished loading IMGT DB " + dbversion)
    if verbose:
        logging.info("Finished loading ALL DB versions")
    gfe_df = pd.DataFrame(gfe_e, columns=":START_ID(ALLELE),:END_ID(ALLELE),imgt_release,:TYPE".split(","))
    seq_df = pd.DataFrame(seq_e, columns=":START_ID(ALLELE),:END_ID(SEQUENCE),imgt_release,accession,:TYPE".split(","))
    seqn_df = pd.DataFrame(seq_n, columns="sequenceId:ID(SEQUENCE),sequence,name,feature:LABEL,rank,length,seq:string[]".split(","))
    allele_df = pd.DataFrame(allele_n, columns="alleleId:ID(ALLELE),name,alleletype:LABEL,locus".split(","))
    group_df = pd.DataFrame(grp_e, columns=":START_ID(ALLELE),:END_ID(ALLELE),imgtdb,:TYPE".split(","))
    cdsn_df = pd.DataFrame(cds_n, columns="cdsId:ID(CDS),name,cdstype:LABEL,cds,protein".split(","))
    trs_df = pd.DataFrame(trs_e, columns=":START_ID(SEQUENCE),:END_ID(CDS),:TYPE".split(","))

    if verbose:
        gfe_es = str(len(gfe_df))
        seq_es = str(len(seq_df))
        seq_ns = str(len(seqn_df))
        all_ns = str(len(allele_df))
        grp_es = str(len(group_df))
        cds_ns = str(len(cdsn_df))
        cds_es = str(len(trs_df))
        logging.info("GFE Edges    = " + gfe_es)
        logging.info("Seq Edges    = " + seq_es)
        logging.info("Group Edges  = " + grp_es)
        logging.info("CDS Edges    = " + cds_es)
        logging.info("Seq Nodes    = " + seq_ns)
        logging.info("CDS Nodes    = " + cds_ns)
        logging.info("Allele Nodes = " + all_ns)

    gfe_df.to_csv(outdir + "/gfe_edges.csv", header=True, index=False)
    seq_df.to_csv(outdir + "/seq_edges.csv", header=True, index=False)
    seqn_df.to_csv(outdir + "/sequence_nodes.csv", header=True, index=False)
    allele_df.to_csv(outdir + "/allele_nodes.csv", header=True, index=False)
    cdsn_df.to_csv(outdir + "/cds_nodes.csv", header=True, index=False)
    group_df.to_csv(outdir + "/group_edges.csv", header=True, index=False)
    trs_df.to_csv(outdir + "/cds_edges.csv", header=True, index=False)

    if verbose:
        logging.info("** Finshed build **")

if __name__ == '__main__':
    """The following will be run if file is executed directly,
    but not if imported as a module"""
    main()


