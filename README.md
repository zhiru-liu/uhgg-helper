# UHGG helper 

This module provides helper codes for dealing with UHGG SNV catalogs.

## Table of Contents
- [Overall workflow](#overall-workflow)
- [Haplotype count table](#haplotype-count-table)
- [Necessary data from UHGG](#necessary-data-from-uhgg)
- [Detailed workflow & key scripts](#detailed-workflow--key-scripts)
    - [Splitting SNV catalog](#splitting-snv-catalog)
    - [Computing core genes](#computing-core-genes)
    - [Compute haplotype count table for pairs of sites](#compute-haplotype-count-table-for-pairs-of-sites)


## Overall workflow
For this project, the goal is to take a given species, specified by the MGNIFY accession id (e.g. `MGYG-HGUT-02492`), and generate a [haplotype count table](#haplotype-count-table) for pairwise-linkage analysis. This is done in the following steps:

1. Identify the accession of the desired species & reference genome in species metadata
2. Download the reference genome & prepare annotation (gene locations & SNV type)
3. Download the SNV catalog of the desired species onto cluster or laptop
4. Split the whole table into tables of individuals genes (to save memory)
5. Compute core genes for the species
6. (Optional) Filter strains to keep only the major clade
7. Process gene tables into haplotype count table (see below) 
    1. Filter strains to keep only “non-redundant” ones
    2. Filter genes to keep only CDS (coding sequence) and “core genes”

## Haplotype count table
Below is an example haplotype count table:

| n11 | n10 | n01 | n00  | ell | type |
|-----|-----|-----|------|-----|------|
|  0  |  2  |  0  | 4692 |  5  |  1   |
|  0  |  2  |  0  | 4691 |  6  |  0   |
|  0  |  2  |  0  | 4691 |  7  |  1   |
|  0  |  2  |  1  | 4690 |  8  |  2   |

Each row is computed for a pair of sites. For simplicity, only bi-allelic sites are considered. In order to annotate the SNV type, only coding regions of non-overlapping genes are considered. The columns are as follows:

- `n11`: number of haplotypes with both SNVs at the two sites (i.e. $n_{AB}$)
- `n10`: number of haplotypes with SNV at the first site but not the second (i.e. $n_{Ab}$)
- `n01` and `n00`: similarly defined. Notably `n00` is the number of haplotypes of with the reference allele at both sites, thus showing the largest number of counts.
- `ell`: distance between two sites (minimum 1 for neighboring sites)
- `type`: 0: syn & syn pair; 1: non-syn & syn pair; 2: non-syn & non-syn pair

Site pairs with similar distances can then be grouped together for linkage analysis. See the estimator module for more details.

## Necessary data from UHGG
We used v1.0 of the UHGG database: http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/

The necessary data are:
- SNV catalogs: found in snv_catalogue/ directory (e.g. `MGYG-HGUT-02492.tsv.tar.lz4`)
- Gene annotations: found in gene_catalogue/ directory (e.g. `MGYG-HGUT-02492.gff`)
- Metadata of non-redundant genomes: `genomes-nr_metadata.tsv`
- Metadata of UHGG species: Supp Table 2 of [Almeida et al. 2021](https://doi.org/10.1038/s41587-020-0603-3); provided in `uhgg_helper/dat/41587_2020_603_MOESM3_ESM.csv`

## Detailed workflow & key scripts 

### Splitting SNV catalog
- Use script `./scripts/split_snvs.py`
- Need the gff file of the species, and the full snv table in `config.SNV_DIR`
- Output to `config.SNV_TABLE_DIR`
- Can annotate the snv table and add two extra columns `snv type` and `degeneracy` using functions in `annotation_utils.py`. The result is stored in `config.ANNOTATED_SNV_DIR`. 

### Computing core genes
- Use script `./scripts/compute_species_core_genes.py`
- Need splitted snv tables, non-redundant genome metadata `config.nr_genome_path`
- Output to `config.CORE_GENE_DIR`

### Compute haplotype count table for pairs of sites
- Use script `./scripts/compute_site_pairs.py`
- Need all information from the previous steps
- Output to `config.SITEPAIR_DIR`

<!-- ### Resulting folder structure
```
├── uhgg
    ├── core_genes
    │   └── MGYG-HGUT-00001
│   ├── genes
│   │   └── MGYG-HGUT-02492.gff
│   ├── reference_genomes
│   │   └── MGYG-HGUT-02492.fna
│   └── snv_tables
│       └── MGYG-HGUT-02492

``` -->

