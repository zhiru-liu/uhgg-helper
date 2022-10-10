import sys
import config
import os
import pandas as pd
import species_diversity_utils, UHGG_utils

if len(sys.argv) !=2:
    print("Usage: python compute_pairwise_distances.py MGNIFY_ACCESSION")
    quit()
accession = sys.argv[1]

print("Processing {}".format(accession))
gff_path = os.path.join(config.GFF_DIR, '{}.gff'.format(accession))
gene_df = UHGG_utils.gff_to_df(gff_path)
grouped_snvs_base = os.path.join(config.SNV_TABLE_DIR, accession)

gene_files = os.listdir(grouped_snvs_base)
header = UHGG_utils.get_SNVs_table_header_from_accession(accession)
genome_names = UHGG_utils.get_genome_names_from_table_header(header)
genomes_metadata = pd.read_csv(config.nr_genome_path, delimiter='\t')
genome_mask = UHGG_utils.get_non_redundant_genome_mask(genomes_metadata, accession, genome_names)

species_diversity_utils.compute_pairwise_distances(accession, grouped_snvs_base, genome_mask, gene_df, header=header)
