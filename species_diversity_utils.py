import numpy as np
import pandas as pd
import os
import UHGG_utils
import config
from datetime import datetime


def compute_pairwise_distances(accession, grouped_snvs_base, genome_mask, gene_df, header=None, debug=False):
    species_pd_dir = os.path.join(config.PD_DIR, '{}'.format(accession))
    os.makedirs(species_pd_dir, exist_ok=True)

    num_genomes = genome_mask.sum()
    cumu_pair_snvs = np.zeros((num_genomes, num_genomes))
    cumu_pair_passed = np.zeros((num_genomes, num_genomes))
    cumu_diff_genes = np.zeros((num_genomes, num_genomes))
    cumu_passed_genes = np.zeros((num_genomes, num_genomes))

    species_core_genes = UHGG_utils.load_core_genes(accession)
    if header is None:
        header = UHGG_utils.get_SNVs_table_header_from_accession(accession)

    processed = 0
    total_sites = 0

    gene_files = os.listdir(grouped_snvs_base)
    print("Total {} genes to process".format(len(gene_files)))
    for gene_file in gene_files:
        gene_file_path = os.path.join(grouped_snvs_base, gene_file)
        items = gene_file.split('.')[0].split('-')
        gene_id = int(items[-1])  # this gene id is assigned by us according to the reference genome's annotation order
        contig = items[0]
        if gene_df.loc[gene_id, 'Type'] != 'CDS':
            continue
        elif gene_id not in species_core_genes:
            # check if gene is a "core gene"
            continue

        # load SNV table
        dat = pd.read_csv(gene_file_path, delimiter='\t', header=None)
        dat.columns = header.split('\t')
        dat.sort_values('Pos', inplace=True)

        # filtering sites with more than two alleles
        biallelic_dat = dat.groupby('Pos').filter(lambda x: x.shape[0] == 1)

        # converting to numpy array for faster arithmetics
        snvs = biallelic_dat.iloc[:, 4:].astype(int).to_numpy()  # use non-annotated snv tables to make sure 4 is correct!
        snvs = snvs[:, genome_mask]

        zeros = (snvs == 0).astype(int)
        ones = (snvs == 1).astype(int)
        covered = (snvs != 255).astype(int)

        pair_snvs = np.dot(zeros.T, ones)  # using dot product to count the number of sites that are different
        pair_snvs = pair_snvs + pair_snvs.T
        pair_passed = np.dot(covered.T, covered)

        total_sites += snvs.shape[0]
        cumu_pair_snvs += pair_snvs
        cumu_pair_passed += pair_passed  # record the number of covered sites between a pair
        cumu_diff_genes += pair_snvs > 0  # record if a pair has non identical genes
        cumu_passed_genes += pair_passed > 0  # record if a gene has any covered sites between a pair

        if processed % 100 == 0:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print("Finished {} at {}".format(processed, current_time))
            if debug and (processed > 100):
                print("Debug mode: exitting")
                break
        processed += 1
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Finished all at {}".format(current_time))

    np.save(os.path.join(species_pd_dir, 'pair_snvs'), cumu_pair_snvs)
    np.save(os.path.join(species_pd_dir, 'pair_passed'), cumu_pair_passed)
    np.save(os.path.join(species_pd_dir, 'diff_genes'), cumu_diff_genes)
    np.save(os.path.join(species_pd_dir, 'passed_genes'), cumu_passed_genes)
