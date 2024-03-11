"""
Compute core genes for a species based on coverage across samples
Core is defined to be genes that are >50% covered in >90% of samples
"""

import pandas as pd
import numpy as np
import sys
import config
import os
import json
import UHGG_utils
import config

if len(sys.argv) !=2:
    print("Usage: python compute_species_core_genes.py MGNIFY_ACCESSION")
    quit()
accession = sys.argv[1]

print("Processing {}".format(accession))
grouped_snvs_base = os.path.join(config.SNV_TABLE_DIR, accession)
genomes_metadata = pd.read_csv(config.nr_genome_path, delimiter='\t')
header = UHGG_utils.get_SNVs_table_header(os.path.join(config.SNV_DIR, '{}_snvs.tsv'.format(accession)))
genome_names = UHGG_utils.get_genome_names_from_table_header(header)

genome_mask = UHGG_utils.get_non_redundant_genome_mask(genomes_metadata, accession, genome_names)
if not os.path.exists(grouped_snvs_base):
    raise RuntimeError("SNV path not found")
coverage_array, gene_files = UHGG_utils.compute_gene_SNV_coverage(grouped_snvs_base, genome_mask, debug=False)
species_core_genes = UHGG_utils.get_species_core_genes(coverage_array, gene_files)

core_gene_root = os.path.join(config.CORE_GENE_DIR, format(accession))
if not os.path.exists(core_gene_root):
    os.makedirs(core_gene_root)
np.save(os.path.join(core_gene_root, 'SNV_coverage'), coverage_array)
json.dump(list(gene_files), open(os.path.join(core_gene_root, 'gene_files.json'.format(accession)), 'w'))
json.dump(list(species_core_genes), open(os.path.join(core_gene_root, 'core_genes.json'), 'w'))
