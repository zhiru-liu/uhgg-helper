"""
This script is used to compute site pairs within and between genes for a given MGNIFY accession.
It takes a single command line argument, which is the MGNIFY accession.
The script processes the SNV files associated with the genes in the accession and generates output files
containing information about the site pairs.

The script performs the following steps:
1. Reads the MGNIFY accession from the command line arguments.
2. Loads necessary configurations and utilities.
3. Creates an output directory for the site pair files.
4. Processes SNV files of each gene in the species.
5. Processes SNV files between randomly selected pairs of genes in the accession.
6. Generates output files containing the haplotype counts between the site pairs.

Note: This script requires the presence of certain directories and files as specified in the 'config' module.
"""

# Rest of the code...
from datetime import datetime
import os
import sys
import config
import UHGG_utils
import pandas as pd
import numpy as np

DEBUG = False
if len(sys.argv) !=2:
    print("Usage: python compute_site_pairs.py MGNIFY_ACCESSION")
    quit()
accession = sys.argv[1]

print("Processing {}".format(accession))
grouped_snvs_base = os.path.join(config.ANNOTATED_SNV_DIR, accession)
gene_files = os.listdir(grouped_snvs_base)

def filename_to_gene_id(filename):
    return int(filename.split('.')[0].split('-')[1])

output_path = os.path.join(config.SITEPAIR_DIR, '{}'.format(accession))
if not os.path.exists(output_path):
    os.makedirs(output_path)
output_path = os.path.join(output_path, '{}.txt'.format(accession))

# write header
# n11: number of haplotypes with both SNVs (i.e. AB)
# n10: number of haplotypes with only the first SNV (i.e. Ab), etc.
# ell: distance between the two sites
# type: 0 for two synonymous SNVs, 1 for one synonymous and one non-synonymous SNV, 2 for two non-synonymous SNVs
with open(output_path, 'w') as f:
    f.write(' '.join(('n11', 'n10', 'n01', 'n00', 'ell', 'type')) + '\n')

processed = 0

# load lots of species information
core_genes = UHGG_utils.load_core_genes(accession)
gff_path = os.path.join(config.GFF_DIR, '{}.gff'.format(accession))
gene_df = UHGG_utils.gff_to_df(gff_path)
print("Total {} genes, {} core_genes".format((gene_df['Type']=='CDS').sum(), len(core_genes)))
snvs_path = os.path.join(config.SNV_DIR, '{}_snvs.tsv'.format(accession))
header = UHGG_utils.get_SNVs_table_header(snvs_path)
genome_names = UHGG_utils.get_genome_names_from_table_header(header)
genomes_metadata = pd.read_csv(config.nr_genome_path, delimiter='\t')
genome_mask = UHGG_utils.get_non_redundant_genome_mask(genomes_metadata, accession, genome_names)

print("Start processing within gene")
for gene_file in gene_files:
    gene_file_path = os.path.join(grouped_snvs_base, gene_file)
    gene_id = filename_to_gene_id(gene_file)
    if processed % 100 == 0:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Finished {} at {}".format(processed, current_time))
    processed += 1
    if gene_df.loc[gene_id, 'Type'] != 'CDS':
        continue
    elif gene_id not in core_genes:
        continue
    UHGG_utils.process_SNV_table_single_gene(gene_file_path, genome_mask, output_path)
    # TODO: eventually add codes here to take care of pairs of SNVs within 4 consecutive genes
    # dealing with contig and gene id is a little bit annoying
    if DEBUG:
        break
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Finished all within gene at {}".format(processed, current_time))

print("Start processing between gene")
n_pairs = 1000
for i in range(n_pairs):
    gene1, gene2 = UHGG_utils.sample_random_pair_of_genes(gene_files, core_genes, gene_df)
    gene_file_path1 = os.path.join(grouped_snvs_base, gene1)
    gene_file_path2 = os.path.join(grouped_snvs_base, gene2)
    if i % 100 == 0:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Finished {} at {}".format(i, current_time))
    UHGG_utils.process_SNV_table_between_genes(gene_file_path1, gene_file_path2, genome_mask, output_path)
    if DEBUG:
        break

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Finished all between gene at {}".format(current_time))

