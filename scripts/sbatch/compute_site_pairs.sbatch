#!/usr/bin/bash
#SBATCH --job-name=site_pairs
#SBATCH --output=out/site_pairs.%j.%a.out
#SBATCH --error=out/site_pairs.%j.%a.err
#SBATCH --time=08:00:00
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --mem=4GB
#SBATCH --array=1-10

module load python/3.6.1
module load py-pandas/1.0.3_py36
export PYTHONPATH=$PYTHONPATH:"${HOME}/uhgg-helper/"

#accession="MGYG-HGUT-02492"
accession=`sed -n "${SLURM_ARRAY_TASK_ID}p" < ../accession_list.txt`
python3 ../compute_site_pairs.py $accession
