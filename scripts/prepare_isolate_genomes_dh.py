import os
import numpy as np
import pandas as pd
from Bio import SeqIO
import argparse
from datetime import datetime
import config
import UHGG_utils, annotation_utils


def get_isolate_genome_mask(genomes_metadata, mgnify_accession, genomes):
    species_rows = genomes_metadata[genomes_metadata['MGnify_accession'] == mgnify_accession]
    species_genomes = species_rows[species_rows['Genome_type']=='Isolate']['Genome']
    species_genomes = np.array(species_genomes)
    return np.isin(genomes, species_genomes)


parser = argparse.ArgumentParser()
parser.add_argument('--accession', type=str, required=True,
                    help='Accession of the species')
parser.add_argument('--debug', action='store_true')
args = parser.parse_args()

accession = args.accession
# accession = 'MGYG-HGUT-02478'
DEBUG = args.debug

# parsing the fasta file for ref seqeunce
input_file = os.path.join(config.REF_GENOME_DIR, "{}.fna".format(accession))
fasta_sequences = SeqIO.parse(open(input_file),'fasta')

all_records = []
for record in fasta_sequences:
    all_records.append(record)
print("In total {} contigs".format(len(all_records)))
sequence = all_records[0].seq

gff_path = os.path.join(config.GFF_DIR, '{}.gff'.format(accession))
gene_df = UHGG_utils.gff_to_df(gff_path)
gene_df.sort_values(by='Start', inplace=True)
cds_genes = gene_df[gene_df['Feature Type']=='CDS']
gene_id_to_row = {row['Gene ID']: row for i, row in gene_df.iterrows()}


# annotate the ref genome
site_pairs = []
for idx, row in cds_genes.iterrows():
    start = row['Start']
    end = row['End']
    assert((end - start + 1)%3 == 0)  # require all coding genes to have the whole reading frames; can improve later
    subseq = sequence[start-1:end]
    types = annotation_utils.annotate_site_types(subseq, str(row['Strand']))
    gene_name = row['Gene ID']
    for i, var_type in enumerate(types):
        site_pairs.append(((gene_name, i+start), var_type))
site_var_dict = dict(site_pairs)

tbl_dir = os.path.join(config.SNV_TABLE_DIR, accession)
snv_tables = os.listdir(tbl_dir)

snvs_path = os.path.join(config.SNV_DIR, '{}_snvs.tsv'.format(accession))
header = UHGG_utils.get_SNVs_table_header(snvs_path)
genome_names = UHGG_utils.get_genome_names_from_table_header(header)
genomes_metadata = pd.read_csv('genomes-nr_metadata.tsv', delimiter='\t')
genome_mask = get_isolate_genome_mask(genomes_metadata, accession, genome_names)

good_genomes = np.array(genome_names)[genome_mask]
contig_name = gene_df.iloc[0, 0]

# initialize all the data arrays
genome_len = len(sequence)
num_isolates = np.sum(genome_mask)
print("{} has {} isolates".format(accession, num_isolates))
snp_array = np.zeros((genome_len, num_isolates))
covered_array = np.zeros((genome_len, num_isolates), dtype=bool)
chromos = np.full(genome_len, contig_name)
locations = np.arange(genome_len) + 1
variants = np.full(genome_len, 'NA')
gene_names = np.full(genome_len, -1, dtype=int)
pvalues = np.zeros(genome_len)

# fill the variant and gene arrays
for gene_id, loc in site_var_dict:
    variants[loc-1] = site_var_dict[(gene_id, loc)]
    gene_names[loc-1] = gene_id


print(datetime.now())
num_processed = 0
multi_allele_sites = []
for filename in snv_tables:
    path = os.path.join(tbl_dir, filename)
    gene_id = int(filename.split('.')[0].split('-')[-1])
    gene_info = gene_id_to_row[gene_id]
    gene_start = gene_info['Start']
    gene_end = gene_info['End']

    dat = pd.read_csv(path, delimiter='\t', header=None)
    all_covered = dat.iloc[:, 4:].iloc[:, genome_mask] != 255
    covered_genome_mask = all_covered.mean(axis=0) > 0.5
    covered_array[gene_start - 1:gene_end, :][:,
    covered_genome_mask] = 1  # as long as covered in snv sites, should be covered in the whole gene

    all_snvs = dat.iloc[:, 4:].iloc[:, genome_mask]
    fracs = (all_snvs == 1).sum(axis=1) / (all_snvs != 255).sum(axis=1).astype(float)
    true_snvs = dat[fracs > 0]

    multi_allelic_site_pos = true_snvs.groupby(1).filter(lambda x: x.shape[0] > 1).iloc[:,
                             1].unique() - 1  # minus 1 to shift to 0 indexing
    # record how many sites are multi allelic among the isolates
    multi_allele_sites.append(len(multi_allelic_site_pos))
    # treat these sites as missing data
    covered_array[multi_allelic_site_pos, :] = 0

    biallelic_dat = true_snvs.groupby(1).filter(lambda x: x.shape[0] == 1).copy()
    biallelic_dat.sort_values(by=1, inplace=True)
    snvs = biallelic_dat.iloc[:, 4:].astype(int).to_numpy()
    snvs = snvs[:, genome_mask]
    loc_mask = biallelic_dat.iloc[:, 1].to_numpy() - 1
    snp_array[loc_mask] = snvs == 1
    covered_array[loc_mask] = snvs != 255

    num_processed += 1
    if num_processed % 200 == 0:
        print("Finished {} genes".format(num_processed))
        print(datetime.now())
        if DEBUG:
            break

datadir = os.path.join(config.DH_DIR, accession)
os.makedirs(datadir, exist_ok=True)

np.save(os.path.join(datadir, 'chromosomes'), chromos)
np.save(os.path.join(datadir, 'locations'), locations)
np.save(os.path.join(datadir, 'variants'), variants)
np.save(os.path.join(datadir, 'gene_names'), gene_names)
np.save(os.path.join(datadir, 'pvalues'), pvalues)
np.save(os.path.join(datadir, 'snp_array'), snp_array)
np.save(os.path.join(datadir, 'covered_array'), covered_array)
np.save(os.path.join(datadir, 'good_genomes'), good_genomes)
