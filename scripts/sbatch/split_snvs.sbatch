#!/usr/bin/bash
#SBATCH --job-name=split_snvs
#SBATCH --output=out/split_snvs_arr.%j.%a.out
#SBATCH --error=out/split_snvs_arr.%j.%a.err
#SBATCH --time=02:00:00
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --mem=4GB
#SBATCH --array=1-4

module load python/3.6.1
module load py-pandas/1.0.3_py36
export PYTHONPATH=$PYTHONPATH:"${HOME}/uhgg-helper/"

#accession="MGYG-HGUT-02506"
accession=`sed -n "${SLURM_ARRAY_TASK_ID}p" < ../pathogen_isolates.txt`
python3 ../split_snvs.py "${GROUP_SCRATCH}/uhgg/snvs/${accession}_snvs.tsv"
