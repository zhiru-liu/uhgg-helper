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
contig_to_seq = {rec.id: rec.seq for rec in all_records}

gff_path = os.path.join(config.GFF_DIR, '{}.gff'.format(accession))
gene_df = UHGG_utils.gff_to_df(gff_path)
cds_genes = gene_df[gene_df['Type']=='CDS']
gene_id_to_row = {row['Gene ID']: row for i, row in gene_df.iterrows()}

snvs_path = os.path.join(config.SNV_DIR, '{}_snvs.tsv'.format(accession))
header = UHGG_utils.get_SNVs_table_header(snvs_path)
genome_names = UHGG_utils.get_genome_names_from_table_header(header)
genomes_metadata = pd.read_csv(config.nr_genome_path, delimiter='\t')
genome_mask = get_isolate_genome_mask(genomes_metadata, accession, genome_names)
good_genomes = np.array(genome_names)[genome_mask]
num_isolates = np.sum(genome_mask)
print("{} has {} isolates".format(accession, num_isolates))

tbl_dir = os.path.join(config.SNV_TABLE_DIR, accession)
snv_tables = os.listdir(tbl_dir)
# group snv tables by contig
contig_to_tables = {}
for filename in snv_tables:
    contig = filename.split('.')[0].split('-')[0]
    if contig not in contig_to_tables:
        contig_to_tables[contig] = []
    contig_to_tables[contig].append(filename)
# save and sort contig by id
all_contigs = contig_to_tables.keys()
all_contigs = sorted(all_contigs, key=lambda x: int(x.split('_')[-1]))

print(datetime.now())
num_processed = 0
# defining lists to hold results for separate contigs
full_chromos = []
full_locations = []
full_variants = []
full_gene_names = []
full_pvalues = []
full_snp_array = []
full_covered_array = []
for contig in all_contigs:
    # contigs that have no snv tables will be skipped
    print("Processing {}".format(contig))
    sequence = contig_to_seq[contig]
    contig_cds_genes = cds_genes[cds_genes['Contig']==contig]
    print("Current contig has {} CDS genes".format(contig_cds_genes.shape[0]))
    # annotate the ref genome contig
    site_pairs = []
    for idx, row in contig_cds_genes.iterrows():
        start = row['Start']
        end = row['End']
        assert((end - start + 1)%3 == 0)  # require all coding genes to have the whole reading frames; can improve later
        subseq = sequence[start-1:end]
        types = annotation_utils.annotate_site_types(subseq, str(row['Strand']))
        gene_name = row['Gene ID']
        for i, var_type in enumerate(types):
            site_pairs.append(((gene_name, i+start), var_type))
    site_var_dict = dict(site_pairs)

    contig_name = contig

    # initialize all the data arrays
    genome_len = len(sequence)
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

    contig_snv_tables = contig_to_tables[contig]
    print("Has in total {} genes".format(len(contig_snv_tables)))
    multi_allele_sites = []
    for filename in contig_snv_tables:
        # now iterate over all the genes in this contig
        path = os.path.join(tbl_dir, filename)
        gene_id = int(filename.split('.')[0].split('-')[-1])
        gene_info = gene_id_to_row[gene_id]
        gene_start = gene_info['Start']
        gene_end = gene_info['End']

        dat = pd.read_csv(path, delimiter='\t', header=None)
        all_covered = dat.iloc[:, 4:].iloc[:, genome_mask] != 255
        covered_genome_mask = all_covered.mean(axis=0) > 0.5
        covered_array[gene_start - 1:gene_end, :][:, covered_genome_mask] = 1  # as long as covered in snv sites, should be covered in the whole gene

        all_snvs = dat.iloc[:, 4:].iloc[:, genome_mask]
        fracs = (all_snvs == 1).sum(axis=1) / (all_snvs != 255).sum(axis=1).astype(float)
        true_snvs = dat[fracs > 0]

        multi_allelic_site_pos = true_snvs.groupby(1).filter(lambda x: x.shape[0] > 1).iloc[:, 1].unique() - 1  # minus 1 to shift to 0 indexing
        # record how many sites are multi allelic among the isolates
        multi_allele_sites.append(len(multi_allelic_site_pos))
        # treat these sites as missing data
        covered_array[multi_allelic_site_pos, :] = 0

        biallelic_dat = true_snvs.groupby(1).filter(lambda x: x.shape[0] == 1).copy()
        biallelic_dat.sort_values(by=1, inplace=True)
        snvs = biallelic_dat.iloc[:, 4:].astype(int).values
        snvs = snvs[:, genome_mask]
        loc_mask = biallelic_dat.iloc[:, 1].values - 1
        snp_array[loc_mask] = snvs == 1
        covered_array[loc_mask] = snvs != 255

        num_processed += 1
	if num_processed % 200 == 0:
	    print("Finished {} genes".format(num_processed))
	    print(datetime.now())
	    if DEBUG:
		break
    full_chromos.append(chromos)
    full_locations.append(locations)
    full_variants.append(variants)
    full_gene_names.append(gene_names)
    full_pvalues.append(pvalues)
    full_snp_array.append(snp_array)
    full_covered_array.append(covered_array)
    if DEBUG:
        # stop after finishing one contig
        break

datadir = os.path.join(config.DH_DIR, accession)
if not os.path.exists(datadir):
    os.makedirs(datadir)

chromos = np.concatenate(full_chromos)
locations = np.concatenate(full_locations)
variants = np.concatenate(full_variants)
gene_names = np.concatenate(full_gene_names)
pvalues = np.concatenate(full_pvalues)
snp_array = np.vstack(full_snp_array)
covered_array = np.vstack(full_covered_array)
print(chromos.shape)
print(snp_array.shape)
assert(snp_array.shape[1]==num_isolates)  # sanity check that the stack direction is correct
assert(snp_array.shape[0]==chromos.shape[0])  # sanity check that the stack direction is correct

np.save(os.path.join(datadir, 'chromosomes'), chromos)
np.save(os.path.join(datadir, 'locations'), locations)
np.save(os.path.join(datadir, 'variants'), variants)
np.save(os.path.join(datadir, 'gene_names'), gene_names)
np.save(os.path.join(datadir, 'pvalues'), pvalues)
np.save(os.path.join(datadir, 'snp_array'), snp_array)
np.save(os.path.join(datadir, 'covered_array'), covered_array)
np.save(os.path.join(datadir, 'good_genomes'), good_genomes)
