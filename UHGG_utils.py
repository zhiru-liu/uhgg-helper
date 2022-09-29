import itertools as it
from datetime import datetime
import os
import pandas as pd
import json
import numpy as np
import random

import config


def gff_to_df(gff_path, save_path=None):
    """
    Parse GFF file downloaded from uhgg into dataframe
    GFF file contains predicted gene annotation
    E.g.:
        http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00001/genome/MGYG-HGUT-00001.gff
    :param gff_path:
    :param save_path: If provided, save the pased file into a csv formatted file for future use
    :return:
    """
    with open(gff_path, 'r') as f:
        lines = f.readlines()
    contigs = []
    sources = []
    types = []
    starts = []
    ends = []
    scores = []
    strands = []
    frames = []
    attributes = []
    for line in lines:
        items = line.split('\t')
        contigs.append(items[0])
        sources.append(items[1])
        types.append(items[2])
        starts.append(int(items[3]))
        ends.append(int(items[4]))
        scores.append(items[5])
        strands.append(items[6])
        frames.append(items[7])
        attributes.append(items[8])
    df = pd.DataFrame(list(zip(contigs, sources, types, starts, ends, scores, strands, frames, attributes)), columns=['Contig', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attribute'])
    df['Gene ID'] = np.arange(df.shape[0])  # We assign gene id by order in gff file. This is an arbitrary choice
    df['UHGG ID'] = df['Attribute'].apply(lambda x: x.split(';')[0].split('=')[-1])
    if save_path is not None:
        df.to_csv(save_path)
    return df


def process_SNVs_table(snvs_path, gene_df, save_path):
    """
    Split the SNV table into smaller ones by gene
    :param snvs_path: path to the big SNV table downloaded from uhgg and decompressed into .tsv
    :param gene_df: parsed by gff_to_df
    :param save_path: Path to folder holding all the individual gene files
    :return:
    """
    linecount = 0
    contig_dfs = {contig: df for contig, df in gene_df.groupby('Contig')}  # faster than slicing df every time changing contig
    missing_contigs = []
    with open(snvs_path, 'r') as f:
        head = f.readline()  # will be useful when processing individual genes
        for line in f:
            items = line[:-2].split('\t')
            contig = items[0]
            # zhiru: newly added because some gff file are missing some contigs (ex species 00060)
            if contig not in contig_dfs:
                if contig not in missing_contigs:
                    missing_contigs.append(contig)
                continue
            contig_df = contig_dfs[contig]
            pos = int(items[1])
            row = if_gene(contig_df, pos)  # sites found in multiple genes will be removed!
            if row is not None:
                # Important: here we set how the gene files should be named
                file_name = os.path.join(save_path, "{}-{}.tsv".format(contig, row['Gene ID'].squeeze()))
                with open(file_name, 'a') as g:
                    g.write(line)
            if linecount % 100000 == 0:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                print("Finished {} at {}".format(linecount, current_time))
            linecount += 1
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Finished all at {}".format(current_time))
    print("Missing these contigs: {}".format(missing_contigs))


def compute_gene_SNV_coverage(grouped_snvs_base, genome_mask, debug=False):
    """
    Compute the estimated coverage for each gene for all non-redundant genomes
    :param grouped_snvs_base: Directory to the separated SNV tables per gene
    :param genome_mask: a boolean mask returned by get_non_redundant_genome_mask; needed to filter genomes
    :return: 2d array of shape num_genes x num_genomes; each element is the fraction of SNV sites that are covered
    """
    gene_files = os.listdir(grouped_snvs_base)
    num_genomes = genome_mask.sum()
    processed = 0
    coverage_array = np.zeros((len(gene_files), num_genomes))
    for gene_file in gene_files:
        gene_file_path = grouped_snvs_base + gene_file
        dat = pd.read_csv(gene_file_path, delimiter='\t', header=None)

        # filtering sites with more than two alleles
        biallelic_dat = dat.groupby(1).filter(lambda x: x.shape[0] == 1)  # column 1 is the positions of SNVs

        # converting to numpy array for faster arithmetics
        snvs = biallelic_dat.iloc[:, 4:].astype(int).to_numpy()
        snvs = snvs[:, genome_mask]

        covered = (snvs != 255).astype(int)
        coverage_array[processed, :] = covered.mean(axis=0)
        if processed % 100 == 0:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print("Finished {} at {}".format(processed, current_time))
            if debug:
                break
        processed += 1
    return coverage_array, gene_files


def get_species_core_genes(coverage_array, gene_files):
    core_gene_mask = (coverage_array > 0.5).mean(axis=1) > 0.9
    gene_ids = [x.split('.')[0].split('-')[-1] for x in gene_files]
    gene_ids = np.array(gene_ids)
    return gene_ids[core_gene_mask]


def compute_core_genes(gene_names, coverage_array, coverage_threshold=0.5, prevelance_threshold=0.9):
    prevelances = np.mean(coverage_array >= coverage_threshold, axis=1)
    return gene_names[prevelances >= prevelance_threshold]

def load_core_genes(species_id):
    # return json.load(open('/Volumes/Botein/uhgg/core_genes/{}/core_genes.json'.format(species_id), 'r'))
    path = os.path.join(config.CORE_GENE_DIR, species_id, 'core_genes.json')
    return json.load(open(path, 'r'))

def get_SNVs_table_header(snvs_path):
    with open(snvs_path, 'r') as f:
        head = f.readline()  # will be useful when processing individual genes
    return head.rstrip('\n')


def get_genome_names_from_table_header(header):
    return header.split('\t')[4:]


def if_gene(df, pos):
    """
    Return the row containing information of the gene of a certain position
    it's possible for a position to be involved in multiple genes; we ignore those cases
    :param df: Gene df as parsed by gff_to_df
    :param pos: Position along the genome
    :return: If unique gene is found, return row; else None
    """
    mask = (df['Start'] <= pos) & (df['End'] >= pos)
    passed = np.sum(mask)
    if passed==1:
        return df[mask]
    else:
        return None


def get_non_redundant_genome_mask(genomes_metadata, mgnify_accession, genomes):
    """
    According to metadata, the total # of genomes can be a few times larger than the # of non-redundant genome
    Therefore, an additional step to filter genomes might be necessary
    :param genomes_metadata: loaded from "genomes-nr_metadata.tsv"
    :param mgnify_accession: species ID defined by the UHGG collection
    :param genomes: List of genome names as in the header of SNV tables
    :return: a boolean mask for unique genomes
    """
    speices_genomes = genomes_metadata[genomes_metadata['MGnify_accession'] == mgnify_accession]['Genome'].to_numpy()
    return np.isin(genomes, speices_genomes)


def sample_random_pair_of_genes(files, core_genes, gene_df):
    """
    Given a list of single gene snv table files, sample two that are both coding, core genes and separated far enough
    :param files: list of single gene snv tables
    :param core_genes: list of core gene ids (assigned by us), computed by get_species_core_genes or loaded by load_core_genes`
    :param gene_df: parsed from gff_to_df
    :return: two gene file names
    """
    while True:
        file1, file2 = random.sample(files, 2)
        id1 = int(file1.split('.')[0].split('-')[-1])
        id2 = int(file2.split('.')[0].split('-')[-1])
        if np.abs(id1-id2) < 15:
            continue
        type1 = gene_df.loc[id1, 'Type']
        type2 = gene_df.loc[id2, 'Type']
        if (type1 != 'CDS') or (type2 != 'CDS'):
            continue
        if (id1 in core_genes) and (id2 in core_genes):
            return file1, file2


def load_biallelic_snvs(snv_path, genome_mask):
    """
    Handy function for loading the snv table of a single gene
    :param snvs_path: path to the snv table file
    :param genome_mask: for filtering the columns (genomes) of the table; removing "redundant genomes"
    :return: two SNV tables of identical to the ref or diff from the ref (true_zeros & true_ones); the reference genome locations
    """
    dat = pd.read_csv(snv_path, delimiter='\t', header=None)
    dat.sort_values(by=1, inplace=True)  # column 1 is the position along reference genome

    # filtering sites with more than two alleles
    # TODO: this can be improved in the future to take care non-biallelic sites
    biallelic_dat = dat.groupby(1).filter(lambda x: x.shape[0] == 1)

    # converting to numpy array for faster arithmetics
    snvs = biallelic_dat.iloc[:, 4:].astype(int).to_numpy()
    snvs = snvs[:, genome_mask]

    zeros = (snvs == 0).astype(int)
    ones = (snvs == 1).astype(int)
    ntot = (snvs != 255).astype(int).sum(axis=1)

    # set the reference to be the major allele
    true_zeros = zeros.copy()
    true_ones = ones.copy()
    to_flip = zeros.sum(axis=1) < ntot * 0.5
    true_zeros[to_flip] = ones[to_flip]
    true_ones[to_flip] = zeros[to_flip]
    return true_zeros, true_ones, biallelic_dat.iloc[:, 1].to_numpy()


def load_biallelic_snvs_annotated(snv_path, genome_mask):
    """
    The same as the above, but with annotation there are two more columns
    """
    dat = pd.read_csv(snv_path, delimiter='\t', header=None)
    dat.sort_values(by=1, inplace=True)  # column 1 is the position along reference genome

    # filtering sites with more than two alleles
    # TODO: this can be improved in the future to take care non-biallelic sites
    biallelic_dat = dat.groupby(1).filter(lambda x: x.shape[0] == 1)

    # converting to numpy array for faster arithmetics
    # snvs = biallelic_dat.iloc[:, 4:].astype(int).to_numpy()
    snvs = biallelic_dat.iloc[:, 6:].astype(int).to_numpy()  # TODO: eventually use headers maybe
    snvs = snvs[:, genome_mask]

    zeros = (snvs == 0).astype(int)
    ones = (snvs == 1).astype(int)
    ntot = (snvs != 255).astype(int).sum(axis=1)

    # set the reference to be the major allele
    true_zeros = zeros.copy()
    true_ones = ones.copy()
    to_flip = zeros.sum(axis=1) < ntot * 0.5
    true_zeros[to_flip] = ones[to_flip]
    true_ones[to_flip] = zeros[to_flip]

    ells = biallelic_dat.iloc[:, 1].to_numpy()
    mut_types = biallelic_dat.iloc[:, 4].to_numpy()
    return true_zeros, true_ones, ells, mut_types


def pair_type(mut_type1, mut_type2):
    """
    0: syn & syn pair
    1: non-syn & syn pair
    2: non-syn & non-syn pair
    :param mut_type1:
    :param mut_type2:
    :return:
    """
    return int(mut_type1 != 's') + (mut_type2 != 's')


def process_SNV_table_single_gene(filepath, genome_mask, output=None):
    """
    Compute genotype counts for all pairs of sites in a single gene
    :param filepath: file path to the gene SNV table
    :param genome_mask: computed by get_non_redundant_genome_mask, for filtering columns of the snv table
    :param output: If provided, pairwise results will be written to the file as a line (n11, n10, n01, n00, ell) separated by whitespace
    :return: List of pairwise results, if output is not provided
    """

    # true_zeros, true_ones, locations = load_biallelic_snvs(filepath, genome_mask)
    true_zeros, true_ones, locations, mut_types = load_biallelic_snvs_annotated(filepath, genome_mask)

    n00s = np.dot(true_zeros, true_zeros.T)  # using dot product to count the number of genomes with 00 genotype
    n10s = np.dot(true_ones, true_zeros.T)
    n11s = np.dot(true_ones, true_ones.T)
    # TODO: in the future, we can also calculate pairs of syn/non-syn in a similar fashion here
    if output is not None:
        with open(output, 'a') as f:
            for i, j in it.combinations(range(true_zeros.shape[0]), 2):
                ell = np.abs(locations[i] - locations[j])
                n00 = str(n00s[i, j])
                n10 = str(n10s[i, j])
                n01 = str(n10s[j, i])
                n11 = str(n11s[i, j])
                type = str(pair_type(mut_types[i], mut_types[j]))
                # res = (n11, n10, n01, n00, str(ell), contig, str(gene_id))
                # shortened output to save memory
                res = (n11, n10, n01, n00, str(ell), type)
                # res = (n11, n10, n01, n00, str(ell))
                f.write(' '.join(res)+'\n')
    else:
        print("No output file provided; aborting")
        return


def process_SNV_table_between_genes(filepath1, filepath2, genome_mask, output=None):
    """
    Analogous to process_SNV_table_single_gene, but now considering pairs of sites from different genes
    :param filepath1: filepath to first snv table
    :param filepath2: ^
    :param genome_mask: for filtering columns of the snv tables
    :param output: If provided, results will be appended to the files as a line (n11, n10, n01, n00, ell)
    :return:
    """
    # true_zeros1, true_ones1, ells1 = load_biallelic_snvs(filepath1, genome_mask)
    # true_zeros2, true_ones2, ells2 = load_biallelic_snvs(filepath2, genome_mask)
    true_zeros1, true_ones1, ells1, mut_types1 = load_biallelic_snvs_annotated(filepath1, genome_mask)
    true_zeros2, true_ones2, ells2, mut_types2 = load_biallelic_snvs_annotated(filepath2, genome_mask)

    n00s = np.dot(true_zeros1, true_zeros2.T)  # using dot product to count the number of genomes with 00 genotype
    n10s = np.dot(true_ones1, true_zeros2.T)
    n01s = np.dot(true_zeros1, true_ones2.T)
    n11s = np.dot(true_ones1, true_ones2.T)

    # TODO: in the future, we can also calculate pairs of syn/non-syn in a similar fashion here
    if output is not None:
        with open(output, 'a') as f:
            for i in range(n00s.shape[0]):
                for j in range(n00s.shape[1]):
                    ell = np.abs(ells1[i] - ells2[j])
                    n00 = str(n00s[i, j])
                    n10 = str(n10s[i, j])
                    n01 = str(n01s[i, j])
                    n11 = str(n11s[i, j])
                    type = str(pair_type(mut_types1[i], mut_types2[j]))
                    # res = (n11, n10, n01, n00, str(ell), contig, str(gene_id))
                    # shortened output to save memory
                    res = (n11, n10, n01, n00, str(ell), type)
                    f.write(' '.join(res) + '\n')
    else:
        print("No output file provided; aborting")
        return
