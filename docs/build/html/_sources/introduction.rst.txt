Introduction
================

Introduction to GFE
-------------------

HLA genotyping via next generation sequencing (NGS) poses challenges for
the use of HLA allele names to analyze and discuss sequence
polymorphism. NGS will identify many new synonymous and non-coding HLA
sequence variants. Allele names identify the types of nucleotide
polymorphism that define an allele (non-synonymous, synonymous and
non-coding changes), but do not describe how polymorphism is distributed
among the individual features (the flanking untranslated regions, exons
and introns) of a gene. Further, HLA alleles cannot be named in the
absence of antigen-recognition domain (ARD) encoding exons. Here, a
system for describing HLA polymorphism in terms of HLA gene features
(GFs) is proposed. This system enumerates the unique nucleotide
sequences for each GF in an HLA gene, and records these in a GF
enumeration notation that allows both more granular dissection of
allele-level HLA polymorphism and the discussion and analysis of GFs in
the absence of ARD-encoding exon sequences.

Approach
--------

The unique sequences in each gene feature of a given HLA gene can be
sequentially numbered, and applied to construct a new name for that
allele consisting of one field for each GF, containing the unique number
for that GF nucleotide sequence and delimited by dashes, prefaced with
the allele name followed by a 'w' (for Workshop) to identify the
provisional nature of this notation.

Example: The HLA-A gene has 8 exons. To represent the sequence of a
variant of this genes involves determination of the 17 features: 3'UTR,
8 exons, 7 introns, 5'UTR. The allele named ``HLA-A*01:01:01:01`` has a
corresponding sequence which can be aligned and annotated and the 17 numbers can be joined by hyphens to produce the
string ``HLA-Aw2-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-4``.

References
----------

`Mack SJ., Human Immunology,
2015 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4674356/>`__
