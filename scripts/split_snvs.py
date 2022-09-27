import pandas as pd
import sys
import config
import os
import UHGG_utils

if len(sys.argv) !=2:
    print("Usage: python split_snv.py TSVFILE")
    quit()
snv_file = sys.argv[1]
accession = snv_file.split('/')[-1].split('_')[0]
print("Processing {}".format(accession))
assert(accession.startswith('MGYG'))

gff_path = os.path.join(config.GFF_DIR, '{}.gff'.format(accession))
gene_df = UHGG_utils.gff_to_df(gff_path)

save_path = os.path.join(config.SNV_TABLE_DIR, accession)
print("Saving to " + save_path)
if not os.path.exists(save_path):
    os.makedirs(save_path)
UHGG_utils.process_SNVs_table(snv_file, gene_df, save_path)
