import os

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

imgt_hla = 'https://www.ebi.ac.uk/ipd/imgt/hla/docs/release.html'
imgt_hla_media_url = 'https://media.githubusercontent.com/media/ANHIG/IMGTHLA/'
imgt_hla_raw_url = 'https://raw.githubusercontent.com/ANHIG/IMGTHLA/'

imgt_kir = 'https://www.ebi.ac.uk/ipd/kir/docs/version.html'
kir_url = 'ftp://ftp.ebi.ac.uk/pub/databases/ipd/kir/KIR.dat'

data_dir = os.path.dirname(__file__) + "/../data/"
# dbversion = 