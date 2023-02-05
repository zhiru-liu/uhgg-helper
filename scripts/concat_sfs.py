import os
import pandas as pd

accession = 'MGYG-HGUT-02492'
dfs = []
for idx in range(4):
    # super manual right now; index range is set by hand when partitioning the gene files
    sfs_path = os.path.expandvars("$GROUP_HOME/uhgg/sfs/{}_{}.csv".format(accession, idx))
    dfs.append(pd.read_csv(sfs_path))
full_df = pd.concat(dfs)
full_df = full_df.drop(full_df.columns[0], axis=1)

sfs_path = os.path.expandvars("$GROUP_HOME/uhgg/sfs/{}_full.csv".format(accession))
full_df.to_csv(sfs_path)
