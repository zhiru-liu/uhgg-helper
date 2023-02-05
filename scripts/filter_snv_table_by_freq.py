import pandas as pd
import numpy as np
import sys
import os
import argparse
import config
import UHGG_utils
import datetime

parser = argparse.ArgumentParser()
parser.add_argument('--accession', type=str, required=True,
                    help='Accession of the species')
parser.add_argument('--savepath', type=str, required=True,
                    help='path to the final file')
parser.add_argument('--savename', type=str, required=True)
parser.add_argument('--fmin', type=float, required=True)
parser.add_argument('--fmax', type=float, required=True)
args = parser.parse_args()

savepath = args.save_path
suffix = args.savename
accession = args.accession
fmin = args.fmin
fmax = args.fmax

# grouped_snvs_base = os.path.join(config.ANNOTATED_SNV_DIR, accession)
grouped_snvs_base = os.path.expandvars("$GROUP_SCRATCH/uhgg_backup/uhgg/snv_tables_annotated/{}".format(accession))
gene_files = os.listdir(grouped_snvs_base)

def filename_to_gene_id(filename):
    return int(filename.split('.')[0].split('-')[1])

# standard loading functions
core_genes = UHGG_utils.load_core_genes(accession)
gff_path = os.path.join(config.GFF_DIR, '{}.gff'.format(accession))
gene_df = UHGG_utils.gff_to_df(gff_path)
snv_dir = os.path.expandvars("$GROUP_SCRATCH/uhgg_backup/uhgg/snvs/")
snvs_path = os.path.join(snv_dir, '{}_snvs.tsv'.format(accession))
header = UHGG_utils.get_SNVs_table_header(snvs_path)
genome_names = UHGG_utils.get_genome_names_from_table_header(header)
genomes_metadata = pd.read_csv(config.nr_genome_path, delimiter='\t')
genome_mask = UHGG_utils.get_non_redundant_genome_mask(genomes_metadata, accession, genome_names)

genes_processed = 0

sfs = UHGG_utils.load_sfs(accession)
good_sfs = sfs[(sfs['Freq'] <= fmax) & (sfs['Freq'] >= fmin)]
genome_gene_pairs = list(zip(list(good_sfs['Genomes']), list(good_sfs['Gene ids']))) # for filtering genes
genome_locs = good_sfs.groupby('Genomes')['Pos'].apply(list).to_dict() # contains all the locations of the sites to keep

savepath = os.path.join(savepath, '{}_filtered_snvs_{}.tsv'.format(accession, suffix))
with open(savepath, 'w') as f:
    for filename in gene_files:
        curr_genome = filename.split('-')[0]
        curr_gene_id = filename.split('-')[1]
        # filter gene type
        if gene_df.loc[curr_gene_id, 'Type'] != 'CDS':
            continue
        elif curr_gene_id not in core_genes:
            continue
        # filter presence of snvs in desired freq
        if (curr_genome, curr_gene_id) not in genome_gene_pairs:
            continue
        filepath = os.path.join(grouped_snvs_base, filename)
        good_locs = genome_locs[curr_genome]
        for line in open(filepath, 'r'):
            items = line.split('\t')
            loc = int(items[1])
            if loc in good_locs:
                f.write(line)
        genes_processed += 1
        if genes_processed % 100 == 0:
            print(genes_processed)
            print(datetime.datetime.now())
