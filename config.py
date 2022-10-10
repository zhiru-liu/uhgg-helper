import os

# define a bunch of important paths; all for cluster use
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
GFF_DIR = os.path.expandvars("$GROUP_HOME/uhgg/genes/")
REF_GENOME_DIR = os.path.expandvars("$GROUP_HOME/uhgg/reference_genomes/")
SNV_DIR = os.path.expandvars("$GROUP_SCRATCH/uhgg/snvs/")
SNV_TABLE_DIR = os.path.expandvars("$GROUP_SCRATCH/uhgg/snv_tables/")
ANNOTATED_SNV_DIR = os.path.expandvars("$GROUP_SCRATCH/uhgg/snv_tables_annotated/")
CORE_GENE_DIR = os.path.expandvars("$GROUP_HOME/uhgg/core_genes/")
SITEPAIR_DIR = os.path.expandvars("$GROUP_HOME/uhgg/site_pairs/")
PD_DIR = os.path.expandvars("$GROUP_HOME/uhgg/pairwise_distances/")

catalog_path = os.path.join(ROOT_DIR, 'dat', '41587_2020_603_MOESM3_ESM.csv')
nr_genome_path = os.path.expandvars("$GROUP_HOME/uhgg/genomes-nr_metadata.tsv")
