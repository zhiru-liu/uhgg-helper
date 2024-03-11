# UHGG helper 

This repo provides helper codes for dealing with UHGG SNV catalogs.

For this project, the goal is to take a given species, specified by the MGNIFY accession id (e.g. MGYG-HGUT-02492), and generate a haplotype count table for pairwise-linkage analysis.
<!-- ## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license) -->

### Haplotype count table
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

$10^{10}$

Write some math in latex: $n_{11}$

### Necessary data from UHGG
The necessary data are organized in the following way:
```
├── uhgg/
│   ├── genes
│   │   └── MGYG-HGUT-02492.gff
│   ├── reference_genomes
│   │   └── MGYG-HGUT-02492.fna
│   └── snv_tables
│       └── MGYG-HGUT-02492
```

### Overall workflow

