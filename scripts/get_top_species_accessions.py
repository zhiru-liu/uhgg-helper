import pandas as pd
import sys
import config

if len(sys.argv) ==1:
    topk = 10
else:
    topk = sys.argv[1]
meta_data = pd.read_csv(config.catalog_path)
meta_data.sort_values('Number of genomes (near-complete)', ascending=False, inplace=True)
top_species = meta_data.head(topk)['MGnify accession'].tolist()
for species in top_species:
    print(species)
