import pandas as pd
import numpy as np
import sys
import os
import argparse
import config
import UHGG_utils
import datetime

parser = argparse.ArgumentParser()
parser.add_argument('--idx', type=int, required=True, 
                    help='Which batch of gene files to process')
parser.add_argument('--accession', type=str, required=True, 
                    help='Accession of the species')
parser.add_argument('--savepath', type=str,
                    help='Path to where filenames are stored')
args = parser.parse_args()

savepath = args.savepath
idx = args.idx
accession = args.accession

split_files = os.path.join(savepath, 'genefiles_{}_{}.txt'.format(accession, idx))

with open(split_files, 'r') as f:
    gene_files = f.read().splitlines()

snv_dir = os.path.expandvars("$GROUP_SCRATCH/uhgg_backup/uhgg/snvs/")
snvs_path = os.path.join(snv_dir, '{}_snvs.tsv'.format(accession))
header = UHGG_utils.get_SNVs_table_header(snvs_path)
genome_names = UHGG_utils.get_genome_names_from_table_header(header)
genomes_metadata = pd.read_csv(config.nr_genome_path, delimiter='\t')
genome_mask = UHGG_utils.get_non_redundant_genome_mask(genomes_metadata, accession, genome_names)

genes_processed = 0
all_ids = []
all_pos = []
all_alts = []
all_tots = []
all_freqs = []

snv_table_dir = os.path.expandvars("$GROUP_SCRATCH/uhgg_backup/uhgg/snv_tables_annotated/")
snv_tables_path = os.path.join(snv_table_dir, accession)

for filename in gene_files:
    table_path = os.path.join(snv_tables_path, filename)
    zeros, ones, ells, types = UHGG_utils.load_biallelic_snvs_annotated(table_path, genome_mask)
    alts = ones.sum(axis=1)
    tots = (ones+zeros).sum(axis=1)
    freq = alts / tots
    
    all_ids.append([filename for i in range(zeros.shape[0])])
    all_pos.append(ells)
    all_freqs.append(freq)
    all_alts.append(alts)
    all_tots.append(tots)
    genes_processed += 1
    if genes_processed % 100 == 0:
        print(genes_processed)
        print(datetime.datetime.now())

ids = np.concatenate(all_ids)
pos = np.concatenate(all_pos)
freqs = np.concatenate(all_freqs)
alts = np.concatenate(all_alts)
tots = np.concatenate(all_tots)
dat_df = pd.DataFrame(data=list(zip(ids, pos, freqs, alts, tots)),
                     columns = ['Gene filename', 'Pos', 'Freq', 'Alt', 'Tot'])

sfs_path = os.path.expandvars("$GROUP_HOME/uhgg/sfs/{}_{}.csv".format(accession, idx))
dat_df.to_csv(sfs_path)
